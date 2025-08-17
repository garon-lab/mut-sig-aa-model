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

import argparse
import sys
import os
import re
from pathlib import Path
import zipfile
import subprocess
from typing import Optional, List

HERE = Path(__file__).resolve().parent

def log(msg: str) -> None:
    print(f"[run_sample_pipeline] {msg}")

def extract_zip(zip_path: Path, dest_dir: Path, keep_extracted: bool) -> Path:
    dest_dir.mkdir(parents=True, exist_ok=True)
    # Extract into a fixed subdir for predictability
    extract_root = dest_dir / "_test_data"
    if extract_root.exists() and keep_extracted:
        log(f"Using existing extracted data at {extract_root}")
        return extract_root
    if extract_root.exists():
        # Clear previous
        for p in sorted(extract_root.rglob('*'), reverse=True):
            try:
                p.unlink()
            except IsADirectoryError:
                pass
        for p in sorted(extract_root.glob('*')):
            if p.is_dir():
                try:
                    p.rmdir()
                except OSError:
                    pass
        try:
            extract_root.rmdir()
        except OSError:
            pass
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

def main() -> None:
    ap = argparse.ArgumentParser(description="Run sample pipeline on test.zip")    
    ap.add_argument("--zip", default="test.zip", help="Path to test data zip (default: test.zip)")
    ap.add_argument("--out_dir", default="results", help="Output directory (default: results)")
    ap.add_argument("--jobs", type=int, default=8, help="Parallel jobs where supported (default: 8)")
    ap.add_argument("--keep-extracted", action="store_true", help="Reuse previously extracted data if present")
    ap.add_argument("--dry-run", action="store_true", help="Only print commands, do not execute")    
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
    dna_manifest = find_manifest(extracted_root, ["dna", "manifest"]) or find_manifest(extracted_root, ["mutect"]) or find_manifest(extracted_root, ["vep"])    
    if dna_manifest is None:
        raise FileNotFoundError("Could not locate a DNA manifest (looked for *dna*manifest*.tsv/.txt/.csv)")

    # 3) Run dna_preprocessor.py (recommended minimal usage, plus make-simplified)
    dna_out = out_dir / "dna"
    dna_out.mkdir(parents=True, exist_ok=True)

    cmd_dna = [
        sys.executable, str(HERE / "dna_preprocessor.py"),
        "--folder", str(extracted_root),
        "--manifest", str(dna_manifest),
        "--out_dir", str(dna_out),
        "--make-simplified",
        "--jobs", str(args.jobs),
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

    # 5) Run multiomic_integration.py (recommended usage)
    integ_out = out_dir / "integration"
    integ_out.mkdir(parents=True, exist_ok=True)
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
    main()
