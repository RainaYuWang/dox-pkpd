"""Validate primary-paper definitions, then retrieve the exact completed baseline.

No model is fitted here. The archived delay-full run supplies unchanged data,
growth parameters and predictions for D, dose beta and both delay stages.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
SOURCE_SHA = "f05543e978e425785fdf918751185cd1a733d300"
RUN_ID = "43137925939"
ARCHIVE_SHA256 = "b39c64c22d8ebe1e9fd6bcd07d42c478c2b105e0d7cfdb1df8f8d1c2d5f593d3"
DEFINITIONS = (
    "time_origin_definition", "D_definition", "theta_definition",
    "initial_cell_number_definition", "fitting_interval_definition",
)


def validate_definitions(paper):
    """Require recorded evidence; checking this schema does not verify a paper."""
    missing = []
    if paper.get("doi") != "10.1038/s41598-017-05902-z":
        missing.append("expected DOI")
    if paper.get("verified") is not True:
        missing.append("verified primary definitions")
    for key in ("origin_in_current_data_hours", "dose_start_in_current_data_hours"):
        value = paper.get(key)
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
            missing.append(key)
    for key in DEFINITIONS:
        if not isinstance(paper.get(key), str) or not paper[key].strip():
            missing.append(key)
        evidence = paper.get("source_evidence", {}).get(key, {})
        for field in ("source", "location", "excerpt"):
            if not isinstance(evidence.get(field), str) or not evidence[field].strip():
                missing.append(f"source_evidence.{key}.{field}")
    if missing:
        raise ValueError("Primary-paper definitions are not ready: " + ", ".join(missing))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check-only", action="store_true")
    args = parser.parse_args()
    paper = json.loads((ROOT / "reports/paper_time_definitions.json").read_text())
    validate_definitions(paper)  # Fail before network access or any MATLAB computation.
    documents = paper.get("source_documents", [])
    expected = {"references/papers/41598_2017_Article_5902.pdf",
                "references/papers/41598_2017_5902_MOESM1_ESM.pdf"}
    if {item["file"] for item in documents} != expected or len(documents) != 2:
        raise ValueError("Record both reviewed primary PDF files and their hashes")
    for item in documents:
        if hashlib.sha256((ROOT / item["file"]).read_bytes()).hexdigest() != item["sha256"]:
            raise ValueError("Primary document hash mismatch: " + item["file"])
    if args.check_only:
        print("Definition fields/evidence and both PDF hashes checked; source interpretation is recorded in the audit.")
        return
    subprocess.run(["git", "fetch", "--no-tags", "--depth=1", "origin", SOURCE_SHA],
                   cwd=ROOT, check=True)
    subprocess.run([sys.executable, "ci/collect_results.py", "--sha", SOURCE_SHA,
                    "--mode", "delay-full", "--run-id", RUN_ID, "--output", "baseline-results"],
                   cwd=ROOT, check=True)
    path = ROOT / "baseline-results" / SOURCE_SHA / (RUN_ID + "-1")
    manifest = json.loads((path / "manifest.json").read_text())
    if manifest["archive_sha256"] != ARCHIVE_SHA256:
        raise RuntimeError("Archived baseline does not match the previously verified results package")
    print(f"Baseline ready: {path}; no old models were re-fitted.")


if __name__ == "__main__":
    try:
        main()
    except ValueError as error:
        print(str(error), file=sys.stderr)
        sys.exit(2)
