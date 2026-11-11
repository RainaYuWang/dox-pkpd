"""Publish an Actions run's outputs through Git without changing the source tree."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tempfile
import time
import zipfile


def git(*args, data=None, env=None, check=True):
    return subprocess.run(['git', *args], input=data, stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE, env=env, check=check)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--phase', choices=['running', 'finished'], required=True)
    parser.add_argument('--setup-outcome', default='pending')
    parser.add_argument('--fit-outcome', default='pending')
    args = parser.parse_args()
    sha = os.environ['GITHUB_SHA']
    run_id = os.environ['GITHUB_RUN_ID']
    attempt = os.environ['GITHUB_RUN_ATTEMPT']
    repo = os.environ['GITHUB_REPOSITORY']
    prefix = f'runs/{sha}/{run_id}-{attempt}/'
    ref = 'refs/heads/pkpd-results'
    manifest = dict(source_sha=sha, run_id=run_id, attempt=attempt,
                    run_url=f'https://github.com/{repo}/actions/runs/{run_id}',
                    phase=args.phase, setup_outcome=args.setup_outcome,
                    fit_outcome=args.fit_outcome, updated_unix=time.time())
    files = []
    for directory in ['joint_results', 'ci-logs']:
        for path in sorted(Path(directory).rglob('*')):
            if path.is_file() and not path.is_symlink() and path.suffix in {
                    '.png', '.csv', '.mat', '.log', '.json', '.txt'}:
                files.append(path)
    status_path = Path('joint_results/ci_status.json')
    if status_path.exists():
        manifest['matlab'] = json.loads(status_path.read_text())
    with tempfile.TemporaryDirectory(prefix='pkpd-publish-') as temp:
        temp = Path(temp)
        archive = temp / 'results.zip'
        with zipfile.ZipFile(archive, 'w', zipfile.ZIP_DEFLATED) as out:
            for path in files:
                out.write(path, path.as_posix())
        if archive.stat().st_size > 90 * 1024**2:
            raise RuntimeError('Result archive exceeds safe Git blob size; use Actions artifacts.')
        payload = archive.read_bytes()
        manifest['archive_sha256'] = hashlib.sha256(payload).hexdigest()
        manifest['files'] = [p.as_posix() for p in files]
        blobs = {'results.zip': payload,
                 'manifest.json': (json.dumps(manifest, indent=2) + '\n').encode()}
        git_env = os.environ.copy()
        git_env.update(GIT_INDEX_FILE=str(temp / 'index'),
                       GIT_AUTHOR_NAME='github-actions[bot]',
                       GIT_AUTHOR_EMAIL='41898282+github-actions[bot]@users.noreply.github.com',
                       GIT_COMMITTER_NAME='github-actions[bot]',
                       GIT_COMMITTER_EMAIL='41898282+github-actions[bot]@users.noreply.github.com')
        for retry in range(4):
            exists = git('ls-remote', '--exit-code', 'origin', ref, check=False)
            if exists.returncode == 0:
                git('fetch', '--no-tags', 'origin', ref)
                parent = git('rev-parse', 'FETCH_HEAD').stdout.decode().strip()
                git('read-tree', parent, env=git_env)
            elif exists.returncode == 2:
                parent = None
                git('read-tree', '--empty', env=git_env)
            else:
                raise RuntimeError('Cannot read results branch: ' + exists.stderr.decode())
            for name, content in blobs.items():
                blob = git('hash-object', '-w', '--stdin', data=content).stdout.decode().strip()
                git('update-index', '--add', '--cacheinfo', '100644', blob,
                    prefix + name, env=git_env)
            tree = git('write-tree', env=git_env).stdout.decode().strip()
            parents = ['-p', parent] if parent else []
            commit = git('commit-tree', tree, *parents, '-m',
                         f'PKPD {sha[:12]} run {run_id}-{attempt}: {args.phase}',
                         env=git_env).stdout.decode().strip()
            pushed = git('push', 'origin', f'{commit}:{ref}', check=False)
            if pushed.returncode == 0:
                print(f'Published {prefix} ({args.phase}, {len(files)} files)')
                return
            if retry == 3:
                raise RuntimeError('Result publication failed: ' + pushed.stderr.decode())
            time.sleep(2)


if __name__ == '__main__':
    main()
