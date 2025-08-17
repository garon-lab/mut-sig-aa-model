"""
RUN SAMPLE PIPELINE
----------------------
Extracts a test dataset zip and runs the recommended usage of:
  1) dna_preprocessor.py
  2) make_manifest.py (build-main mode)
  3) multiomic_integration.py

Defaults assume a `test.zip` laid out with subfolders like dna/, rna/, ch3/, cn/, protein/
and GDC-style manifests.

Usage:
    python run_sample_pipeline.py         --zip test.zip         --out_dir results         [--jobs 8] [--keep-extracted]

What it does:
- Extracts the zip to ./_test_data (unless already extracted or --keep-extracted).
- Auto-discovers a DNA manifest matching *dna*manifest*.tsv|.txt in the extracted tree.
- Runs dna_preprocessor.py with the discovered manifest and writes outputs to <out_dir>/dna.
- Writes a case list (<out_dir>/dna/case_ids.txt) via --make-simplified (used by integration).
- Runs make_manifest.py --build-main to create a unified manifest set in <out_dir>/manifests.
- Runs multiomic_integration.py with --folder <extracted_root>, --manifest <case_ids.txt>, --out_dir <out_dir>/integration, --step all.

"""

#!/usr/bin/env python3
import argparse
import shutil
import sys
import os
import re
from pathlib import Path
import zipfile
import subprocess
from typing import Optional, List
import logging
import multiprocessing
from concurrent.futures import ProcessPoolExecutor, as_completed
import tempfile

HERE = Path(__file__).resolve().parent

# ---------- Logging ----------
logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")
LOGGER = logging.getLogger("run_sample_pipeline")

def log(msg: str) -> None:
    print(f"[run_sample_pipeline] {msg}")

# ---------- Utilities ----------
def extract_zip(zip_path: Path, dest_dir: Path, keep_extracted: bool) -> Path:
    dest_dir.mkdir(parents=True, exist_ok=True)
    extract_root = dest_dir / "_test_data"
    if extract_root.exists():
        if keep_extracted:
            log(f"Using existing extracted data at {extract_root}")
            return extract_root
        else:
            shutil.rmtree(extract_root, ignore_errors=True)
    extract_root.mkdir(parents=True, exist_ok=True)
    log(f"Extracting {zip_path} -> {extract_root}")
    with zipfile.ZipFile(zip_path, 'r') as zf:
        zf.extractall(extract_root)
    return extract_root

def find_manifest(root: Path, pattern_keywords: List[str]) -> Optional[Path]:
    # Search for files containing ALL keywords in their name, case-insensitive
    candidates = []
    for p in root.rglob("*.*"):
        name = p.name.lower()
        if any(name.endswith(ext) for ext in (".tsv", ".txt", ".csv")):
            if all(kw in name for kw in pattern_keywords):
                candidates.append(p)
    # Prefer .tsv > .txt > .csv, shortest path first
    ext_rank = {".tsv": 0, ".txt": 1, ".csv": 2}
    candidates.sort(key=lambda p: (ext_rank.get(p.suffix.lower(), 3), len(str(p))))
    return candidates[0] if candidates else None

def run(cmd: List[str], dry_run: bool = False) -> None:
    log("$ " + " ".join(cmd))
    if dry_run:
        return
    result = subprocess.run(cmd, stdout=sys.stdout, stderr=sys.stderr)
    if result.returncode != 0:
        raise SystemExit(result.returncode)

def read_case_ids(case_list_path: Path) -> List[str]:
    if not case_list_path.exists():
        raise FileNotFoundError(f"Case list not found: {case_list_path}")
    ids: List[str] = []
    with case_list_path.open() as f:
        for line in f:
            s = line.strip()
            if s:
                ids.append(s)
    if not ids:
        raise ValueError(f"No case IDs found in {case_list_path}")
    LOGGER.info("Loaded %d case IDs from %s", len(ids), case_list_path)
    return ids

def write_single_id_manifest(dest_dir: Path, case_id: str) -> Path:
    """
    Create a tiny manifest file acceptable to multiomic_integration.py that contains
    exactly one ID. We keep it as plain text with just the ID per line, which matches
    the repo's common 'case list' pattern.
    """
    dest_dir.mkdir(parents=True, exist_ok=True)
    path = dest_dir / f"{case_id}.txt"
    with path.open("w") as f:
        f.write(case_id + "\n")
    return path

def run_integration_for_case(
    case_id: str,
    extracted_root: Path,
    base_out: Path,
    dry_run: bool,
) -> str:
    """
    Per-sample worker: creates a small one-ID manifest and calls multiomic_integration.py
    directing outputs into a per-case subdirectory.
    """
    single_out = base_out / case_id
    single_out.mkdir(parents=True, exist_ok=True)
    tmp_manifest_dir = base_out / "_tmp_manifests"
    tmp_manifest = write_single_id_manifest(tmp_manifest_dir, case_id)

    cmd_integ = [
        sys.executable, str(HERE / "multiomic_integration.py"),
        "--folder", str(extracted_root),
        "--manifest", str(tmp_manifest),
        "--out_dir", str(single_out),
        "--step", "all",
    ]
    run(cmd_integ, dry_run)
    return case_id

# ---------- CLI / Main ----------
def main() -> None:
    ap = argparse.ArgumentParser(description="Run sample pipeline on test.zip (parallel-enabled)")
    ap.add_argument("--zip", default="test.zip", help="Path to test data zip (default: test.zip)")
    ap.add_argument("--out_dir", default="results", help="Output directory (default: results)")
    ap.add_argument("--jobs", type=int, default=max(1, (multiprocessing.cpu_count() or 1) - 1),
                    help="Parallel jobs where supported (default: CPU count minus one)")
    ap.add_argument("--keep-extracted", action="store_true", help="Reuse previously extracted data if present")
    ap.add_argument("--dry-run", action="store_true", help="Only print commands, do not execute")
    ap.add_argument("--per-sample-integration", action="store_true",
                    help="Run multiomic integration per case in parallel using --jobs")
    args = ap.parse_args()

    dry_run = args.dry_run

    zip_path = Path(args.zip).resolve()
    if not zip_path.exists():
        # Try to resolve relative to script dir
        alt = (HERE / args.zip).resolve()
        if alt.exists():
            zip_path = alt
        else:
            raise FileNotFoundError(f"Zip not found: {args.zip}")

    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    # 1) Extract
    extracted_root = extract_zip(zip_path, HERE, keep_extracted=args.keep_extracted)

    # 2) Discover DNA manifest
    dna_manifest = (
        find_manifest(extracted_root, ["dna", "manifest"])
        or find_manifest(extracted_root, ["mutect"])
        or find_manifest(extracted_root, ["vep"])
    )
    if dna_manifest is None:
        raise FileNotFoundError("Could not locate a DNA manifest (looked for *dna*manifest*.tsv/.txt/.csv)")

    # 3) Run dna_preprocessor.py (recommended minimal usage, plus make-simplified)
    dna_out = out_dir / "dna"
    dna_out.mkdir(parents=True, exist_ok=True)

    cmd_dna = [
        sys.argv[0] if sys.executable is None else sys.executable, str(HERE / "dna_preprocessor.py"),
        "--folder", str(extracted_root),
        "--manifest", str(dna_manifest),
        "--out_dir", str(dna_out),
        "--make-simplified",
        #"--jobs", str(args.jobs),
    ]
    run(cmd_dna, dry_run)

    # Case list path (default from README: <out_dir>/case_ids.txt)
    case_list = dna_out / "case_ids.txt"
    if not case_list.exists():
        # Fallback: try common names
        alt = list(dna_out.glob("*case*id*.txt"))
        if alt:
            case_list = alt[0]
        else:
            raise FileNotFoundError("Case list not found after dna_preprocessor. Expected results/dna/case_ids.txt")

    # 4) Run make_manifest.py --build-main
    mani_out = out_dir / "manifests"
    mani_out.mkdir(parents=True, exist_ok=True)
    cmd_manifest = [
        sys.executable, str(HERE / "make_manifest.py"),
        "--build-main",
        "--project_root", str(extracted_root),
        "--main_out", str(mani_out),
    ]
    run(cmd_manifest, dry_run)

    # 5) Run multiomic_integration.py
    integ_out = out_dir / "integration"
    integ_out.mkdir(parents=True, exist_ok=True)

    if args.per_sample_integration:
        # Parallel per-case execution using ProcessPoolExecutor
        ids = read_case_ids(case_list)
        LOGGER.info("Per-sample integration enabled with --jobs=%d", args.jobs)
        failed = 0
        if args.jobs <= 1 or len(ids) <= 1:
            for sid in ids:
                try:
                    run_integration_for_case(sid, extracted_root, integ_out, dry_run)
                except Exception as e:
                    failed += 1
                    LOGGER.exception("Sample %s failed: %s", sid, e)
        else:
            with ProcessPoolExecutor(max_workers=args.jobs) as ex:
                futures = {
                    ex.submit(run_integration_for_case, sid, extracted_root, integ_out, dry_run): sid
                    for sid in ids
                }
                for fut in as_completed(futures):
                    sid = futures[fut]
                    try:
                        fut.result()
                    except Exception as e:
                        failed += 1
                        LOGGER.exception("Sample %s failed: %s", sid, e)
        if failed:
            LOGGER.warning("Integration completed with %d failed sample(s).", failed)
    else:
        # Original behavior: one integration run over all cases
        cmd_integ = [
            sys.executable, str(HERE / "multiomic_integration.py"),
            "--folder", str(extracted_root),
            "--manifest", str(case_list),
            "--out_dir", str(integ_out),
            "--step", "all",
        ]
        run(cmd_integ, dry_run)

    log("Pipeline complete.")
    log(f"Outputs:\n- DNA: {dna_out}\n- Manifests: {mani_out}\n- Integration: {integ_out}")

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        LOGGER.exception("Pipeline failed: %s", e)
        sys.exit(1)

