"""Fetch completed results for an exact source commit through the Git proxy."""
import argparse
import hashlib
import io
import json
from pathlib import Path
import subprocess
import sys
import time
import zipfile


def git(*args, check=True):
    return subprocess.run(['git', *args], stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE, check=check)


def matches_requested_run(item, mode=None, run_id=None):
    return ((mode is None or item.get('mode', item.get('matlab', {}).get('mode')) == mode)
            and (run_id is None or str(item['run_id']) == str(run_id)))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sha', default='HEAD')
    parser.add_argument('--mode', help='Require this run mode; prevents quick/full result confusion')
    parser.add_argument('--run-id', help='Optionally select an exact Actions run')
    parser.add_argument('--output', default='cloud-results')
    parser.add_argument('--wait', type=int, default=0, help='Maximum seconds to poll; 0 checks once')
    args = parser.parse_args()
    sha = git('rev-parse', '--verify', args.sha + '^{commit}').stdout.decode().strip()
    deadline = time.monotonic() + args.wait
    while True:
        found = git('ls-remote', '--exit-code', 'origin', 'refs/heads/pkpd-results', check=False)
        if found.returncode not in (0, 2):
            raise RuntimeError(found.stderr.decode())
        if found.returncode == 0:
            git('fetch', '--quiet', '--no-tags', 'origin', 'refs/heads/pkpd-results')
            tree = git('rev-parse', 'FETCH_HEAD').stdout.decode().strip()
            paths = git('ls-tree', '-r', '--name-only', tree, f'runs/{sha}/').stdout.decode().splitlines()
            manifests = []
            for path in paths:
                if path.endswith('/manifest.json'):
                    item = json.loads(git('show', f'{tree}:{path}').stdout)
                    if item['source_sha'] != sha:
                        raise RuntimeError('Result source SHA mismatch')
                    if matches_requested_run(item, args.mode, args.run_id):
                        manifests.append((int(item['run_id']), int(item['attempt']), path, item))
            if manifests:
                _, _, path, manifest = max(manifests)
                if manifest['phase'] == 'finished':
                    archive_path = path.removesuffix('manifest.json') + 'results.zip'
                    payload = git('show', f'{tree}:{archive_path}').stdout
                    assert hashlib.sha256(payload).hexdigest() == manifest['archive_sha256']
                    output = Path(args.output) / sha / f"{manifest['run_id']}-{manifest['attempt']}"
                    output.mkdir(parents=True, exist_ok=True)
                    (output / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
                    (output / 'results.zip').write_bytes(payload)
                    with zipfile.ZipFile(io.BytesIO(payload)) as z:
                        for name in z.namelist():
                            if Path(name).is_absolute() or '..' in Path(name).parts:
                                raise RuntimeError('Unsafe path in result archive')
                        z.extractall(output)
                    print(json.dumps(manifest, indent=2), flush=True)
                    print(f'Downloaded and verified results: {output}', flush=True)
                    return 0 if (manifest['fit_outcome'] == 'success' and
                                 manifest.get('matlab', {}).get('status') == 'success') else 1
                print(f"Running: {manifest['run_url']}", flush=True)
            else:
                print(f'No result record for {sha[:12]} yet (possibly queued).', flush=True)
        else:
            print('Results branch not yet created (runner may be queued).', flush=True)
        if time.monotonic() >= deadline:
            return 2
        time.sleep(min(30, max(0, deadline - time.monotonic())))


if __name__ == '__main__':
    sys.exit(main())
