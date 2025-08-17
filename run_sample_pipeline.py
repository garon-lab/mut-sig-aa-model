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
import sys
from pathlib import Path
import zipfile
import subprocess
import shutil
import multiprocessing
from concurrent.futures import ProcessPoolExecutor, as_completed

HERE = Path(__file__).resolve().parent

# ----- Logging -----
def log(msg: str) -> None:
    print(f"[run_sample_pipeline] {msg}")

# ----- Utilities -----
def run(cmd, dry_run=False):
    log("$ " + " ".join(map(str, cmd)))
    if dry_run:
        return
    result = subprocess.run(cmd, stdout=sys.stdout, stderr=sys.stderr)
    if result.returncode != 0:
        raise SystemExit(result.returncode)

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

def read_case_ids(case_list_path: Path):
    case_list_path = case_list_path.expanduser().resolve()
    if not case_list_path.exists():
        raise FileNotFoundError(f"Case list not found: {case_list_path}")
    ids = [ln.strip() for ln in case_list_path.read_text().splitlines() if ln.strip()]
    if not ids:
        raise ValueError(f"No case IDs in {case_list_path}")
    return ids

def write_single_id_manifest(dest_dir: Path, case_id: str) -> Path:
    dest_dir.mkdir(parents=True, exist_ok=True)
    p = dest_dir / f"{case_id}.txt"
    p.write_text(case_id + "\n")
    return p

def run_integration_for_case(case_id: str, extracted_root: Path, base_out: Path, dry_run: bool):
    single_out = base_out / case_id
    single_out.mkdir(parents=True, exist_ok=True)
    tmp_manifest = write_single_id_manifest(base_out / "_tmp_manifests", case_id)
    cmd_integ = [
        sys.executable, str(HERE / "multiomic_integration.py"),
        "--folder", str(extracted_root),
        "--manifest", str(tmp_manifest),
        "--out_dir", str(single_out),
        "--step", "all",
    ]
    run(cmd_integ, dry_run)
    return case_id

# ----- CLI / Main -----
def main():
    ap = argparse.ArgumentParser(description="Run sample pipeline on test.zip")
    ap.add_argument("--zip", default="test.zip", help="Path to test data zip (default: test.zip)")
    ap.add_argument("--out_dir", default="results", help="Output directory (default: results)")
    ap.add_argument("--jobs", type=int,
                    default=max(1, (multiprocessing.cpu_count() or 1) - 1),
                    help="Parallel jobs where supported (default: CPU count minus one)")
    ap.add_argument("--keep-extracted", action="store_true",
                    help="Reuse previously extracted data if present")
    ap.add_argument("--dry-run", action="store_true",
                    help="Only print commands, do not execute")
    ap.add_argument("--per-sample-integration", action="store_true",
                    help="Run multiomic integration per case in parallel using --jobs")
    args = ap.parse_args()

    dry_run = args.dry_run

    zip_path = Path(args.zip).expanduser().resolve()
    if not zip_path.exists():
        alt = (HERE / args.zip).expanduser().resolve()
        if alt.exists():
            zip_path = alt
        else:
            raise FileNotFoundError(f"Zip not found: {args.zip}")

    out_dir = Path(args.out_dir).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    # 1) Extract test data
    extracted_root = extract_zip(zip_path, HERE, keep_extracted=args.keep_extracted)

    # 2) dna_preprocessor.py
    dna_out = out_dir / "dna"
    dna_out.mkdir(parents=True, exist_ok=True)

    dna_manifest = extracted_root / "dna_manifest.tsv"
    if not dna_manifest.exists():
        candidates = []
        for pat in ("*dna*manifest*.tsv", "*dna*manifest*.txt",
                    "*mutect*.tsv", "*mutect*.txt",
                    "*vep*.tsv", "*vep*.txt"):
            candidates.extend(extracted_root.rglob(pat))
        if candidates:
            candidates.sort(key=lambda p: (len(p.suffix), len(str(p))))
            dna_manifest = candidates[0]
        else:
            raise FileNotFoundError("Could not locate a DNA manifest inside the extracted test data.")

    cmd_dna = [
        sys.executable, str(HERE / "dna_preprocessor.py"),
        "--folder", str(extracted_root),
        "--manifest", str(dna_manifest),
        "--out_dir", str(dna_out),
        "--make-simplified",
        "--jobs", str(args.jobs),
    ]
    run(cmd_dna, dry_run)

    # 3) make_manifest.py
    mani_out_dir = out_dir / "manifests"
    mani_out_dir.mkdir(parents=True, exist_ok=True)
    main_manifest = mani_out_dir / "main_manifest.tsv"

    cmd_manifest = [
        sys.executable, str(HERE / "make_manifest.py"),
        "--build-main",
        "--project_root", str(extracted_root),
        "--main_out", str(main_manifest),
    ]
    run(cmd_manifest, dry_run)

    # 4) multiomic_integration.py
    integ_out = out_dir / "integration"
    integ_out.mkdir(parents=True, exist_ok=True)

    case_list = dna_out / "case_ids.txt"
    if not case_list.exists():
        alts = list(dna_out.glob("*case*id*.txt"))
        if alts:
            case_list = alts[0]
        else:
            raise FileNotFoundError("Case list not found after dna_preprocessor (expected results/dna/case_ids.txt).")

    if args.per_sample_integration:
        ids = read_case_ids(case_list)
        log(f"Per-sample integration enabled with --jobs={args.jobs} for N={len(ids)} cases")
        failed = 0
        if args.jobs <= 1 or len(ids) <= 1 or dry_run:
            for sid in ids:
                try:
                    run_integration_for_case(sid, extracted_root, integ_out, dry_run)
                except Exception as e:
                    failed += 1
                    log(f"[WARN] Sample {sid} failed: {e}")
        else:
            with ProcessPoolExecutor(max_workers=args.jobs) as ex:
                futures = {ex.submit(run_integration_for_case, sid, extracted_root, integ_out, dry_run): sid for sid in ids}
                for fut in as_completed(futures):
                    sid = futures[fut]
                    try:
                        fut.result()
                    except Exception as e:
                        failed += 1
                        log(f"[WARN] Sample {sid} failed: {e}")
        if failed:
            log(f"[WARN] Integration completed with {failed} failed sample(s).")
    else:
        cmd_integ = [
            sys.executable, str(HERE / "multiomic_integration.py"),
            "--folder", str(extracted_root),
            "--manifest", str(case_list),
            "--out_dir", str(integ_out),
            "--step", "all",
        ]
        run(cmd_integ, dry_run)

    log("Pipeline complete.")
    log(f"Outputs:\n- DNA: {dna_out}\n- Manifests: {mani_out_dir}\n- Integration: {integ_out}")

if __name__ == "__main__":
    main()
