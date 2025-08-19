#!/usr/bin/env python3
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
    """Extract zip into <dest_dir>/_test_data and return that path."""
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
    ids = [ln.strip() for ln in case_list_path.read_text().splitlines()
           if ln.strip() and not ln.strip().startswith("#")]
    if not ids:
        raise ValueError(f"No case IDs in {case_list_path}")
    return ids

def write_case_ids(case_list_path: Path, ids):
    case_list_path.parent.mkdir(parents=True, exist_ok=True)
    case_list_path.write_text("".join(f"{i}\n" for i in ids))

def write_single_id_manifest(dest_dir: Path, case_id: str) -> Path:
    dest_dir.mkdir(parents=True, exist_ok=True)
    p = dest_dir / f"{case_id}.txt"
    p.write_text(case_id + "\n")
    return p

def _assert_reference_here(reference_zip: Path) -> None:
    if not reference_zip.exists():
        raise FileNotFoundError(
            f"reference.zip not found at {reference_zip}. "
            "Place your reference.zip next to this script or pass a different path."
        )

def _autodiscover_ch3_map(ch3_arg: str) -> Path | None:
    """
    If user provided --ch3_map and it exists, use it.
    Otherwise return None so multiomic_integration.py will rely on --ref_zip (reference.zip).
    """
    if ch3_arg:
        p = Path(ch3_arg).expanduser().resolve()
        if p.exists():
            log(f"Using external CH3 map: {p}")
            return p
        log(f"[WARN] CH3 map not found at {p}; will rely on reference.zip instead")
    else:
        log("No external CH3 map provided; will rely on reference.zip")
    return None

def run_integration_for_case(case_id: str, extracted_root: Path, base_out: Path,
                             ch3_map_path: Path | None, dry_run: bool):
    """Per-sample integration helper (used when --per-sample-integration is on)."""
    single_out = base_out / case_id
    single_out.mkdir(parents=True, exist_ok=True)
    tmp_manifest = write_single_id_manifest(base_out / "_tmp_manifests", case_id)

    # reference + manifests based on outputs
    reference_zip = (HERE / "reference.zip")
    _assert_reference_here(reference_zip)

    mani_out_dir = base_out.parent / "manifests"

    # Raw modality roots under extracted test set
    raw_rna = extracted_root / "test" / "rna"
    raw_ch3 = extracted_root / "test" / "ch3"
    raw_cn  = extracted_root / "test" / "cn"
    raw_pro = extracted_root / "test" / "protein"

    # DNA outputs produced by dna_preprocessor
    dna_out_dir = base_out.parent / "dna" / "dna"

    # Manifests written by make_manifest.py
    rna_manifest = mani_out_dir / "rna_manifest.tsv"
    ch3_manifest = mani_out_dir / "ch3_manifest.tsv"
    cn_manifest  = mani_out_dir / "cn_manifest.tsv"

    cmd_integ = [
        sys.executable, str(HERE / "multiomic_integration.py"),
        "--folder", str(extracted_root),
        "--manifest", str(tmp_manifest),
        "--out_dir", str(single_out),
        "--step", "all",
        "--ref_zip", str(reference_zip),
        "--input_dna_dir", str(dna_out_dir),
        "--input_rna_dir", str(raw_rna),
        "--input_ch3_dir", str(raw_ch3),
        "--input_cn_dir",  str(raw_cn),
        "--input_protein_dir", str(raw_pro),
        "--rna-manifest", str(rna_manifest),
        "--ch3-manifest", str(ch3_manifest),
        "--cn-manifest",  str(cn_manifest),
        "--ensg_join_mode", "core",
    ]
    if ch3_map_path:
        cmd_integ += [
            "--ch3_map", str(ch3_map_path),
            "--ch3_probe_col", "IllmnID",
            "--ch3_symbol_col", "UCSC_RefGene_Name",
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
    ap.add_argument("--ch3_map", default="",
                    help="Path to CH3 probe→gene map (IllmnID,UCSC_RefGene_Name). If omitted, the reference.zip will be used.")
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

    # 1) Extract test data **into out_dir** (not the repo)
    extracted_root = extract_zip(zip_path, out_dir, keep_extracted=args.keep_extracted)

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
            preferred = {"gdc-vep.tsv", "dna_manifest.tsv"}
            candidates.sort(key=lambda p: (0 if p.name in preferred else 1, len(str(p))))
            dna_manifest = candidates[0]
        else:
            raise FileNotFoundError("Could not locate a DNA manifest inside the extracted test data.")

    cmd_dna = [
        sys.executable, str(HERE / "dna_preprocessor.py"),
        "--folder", str(extracted_root),
        "--manifest", str(dna_manifest),
        "--out_dir", str(dna_out),
        "--make-simplified",
    ]
    run(cmd_dna, dry_run)

    # ensure created dna files are available under extracted_root/dna for convenience
    produced = dna_out / "dna"
    dest = extracted_root / "dna"
    dest.mkdir(parents=True, exist_ok=True)
    for f in produced.glob("*.csv"):
        shutil.copy2(f, dest / f.name)

    # If the case list missed any produced DNA CSVs, union them back in (ensures case-01 isn't skipped)
    case_list = dna_out / "case_ids.txt"
    existing_ids = set()
    if case_list.exists():
        try:
            existing_ids = set(read_case_ids(case_list))
        except Exception:
            existing_ids = set()
    produced_ids = {p.stem for p in (dna_out / "dna").glob("case-*.csv")}
    union_ids = sorted(existing_ids | produced_ids)
    if union_ids:
        write_case_ids(case_list, union_ids)
        if produced_ids - existing_ids:
            log(f"Augmented case_ids.txt with: {', '.join(sorted(produced_ids - existing_ids))}")

    # 3) make_manifest.py
    mani_out_dir = out_dir / "manifests"
    mani_out_dir.mkdir(parents=True, exist_ok=True)

    cmd_manifest = [
        sys.executable, str(HERE / "make_manifest.py"),
        "--build-main",
        "--project_root", str(extracted_root),
        "--out_dir", str(mani_out_dir),
    ]
    run(cmd_manifest, dry_run)

    # 4) multiomic_integration.py (batch)
    integ_out = out_dir / "integration"
    integ_out.mkdir(parents=True, exist_ok=True)

    if not case_list.exists():
        raise FileNotFoundError("Case list not found after dna_preprocessor (expected results/dna/case_ids.txt).")

    # Validate reference + (optional) CH3 map
    reference_zip = (HERE / "reference.zip")
    _assert_reference_here(reference_zip)
    ch3_map_path = _autodiscover_ch3_map(args.ch3_map)

    # Build explicit refs and manifests based on outputs and test layout
    raw_rna = extracted_root / "test" / "rna"
    raw_ch3 = extracted_root / "test" / "ch3"
    raw_cn  = extracted_root / "test" / "cn"
    raw_pro = extracted_root / "test" / "protein"
    dna_out_dir = dna_out / "dna"
    rna_manifest = mani_out_dir / "rna_manifest.tsv"
    ch3_manifest = mani_out_dir / "ch3_manifest.tsv"
    cn_manifest  = mani_out_dir / "cn_manifest.tsv"

    if args.per_sample_integration:
        ids = read_case_ids(case_list)
        log(f"Per-sample integration enabled with --jobs={args.jobs} for N={len(ids)} cases")
        failed = 0
        if args.jobs <= 1 or len(ids) <= 1 or dry_run:
            for sid in ids:
                try:
                    run_integration_for_case(sid, extracted_root, integ_out, ch3_map_path, dry_run)
                except Exception as e:
                    failed += 1
                    log(f"[WARN] Sample {sid} failed: {e}")
        else:
            with ProcessPoolExecutor(max_workers=args.jobs) as ex:
                futures = {ex.submit(run_integration_for_case, sid, extracted_root, integ_out, ch3_map_path, dry_run): sid for sid in ids}
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
            "--ref_zip", str(reference_zip),
            "--input_dna_dir", str(dna_out_dir),
            "--input_rna_dir", str(raw_rna),
            "--input_ch3_dir", str(raw_ch3),
            "--input_cn_dir",  str(raw_cn),
            "--input_protein_dir", str(raw_pro),
            "--rna-manifest", str(rna_manifest),
            "--ch3-manifest", str(ch3_manifest),
            "--cn-manifest",  str(cn_manifest),
            "--ensg_join_mode", "core",
        ]
        if ch3_map_path:
            cmd_integ += [
                "--ch3_map", str(ch3_map_path),
                "--ch3_probe_col", "IllmnID",
                "--ch3_symbol_col", "UCSC_RefGene_Name",
            ]
        run(cmd_integ, dry_run)

    log("Pipeline complete.")
    log(f"Outputs:\n- DNA: {dna_out}\n- Manifests: {mani_out_dir}\n- Integration: {integ_out}")

if __name__ == "__main__":
    main()
