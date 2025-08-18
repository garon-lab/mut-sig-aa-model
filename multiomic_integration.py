#!/usr/bin/env python3
"""
MULTIOMIC INTEGRATION

This pipeline integrates multiple omics layers (DNA, RNA, methylation, protein, and copy number)
for each case-id listed in a manifest. Can be used with an existing unifed manifest built with make_manifest.py or separate manifests. Only the final protein files (with copy number added) are retained in the output directory. Note protein files must be preprocessed and in the format {Case-ID}.csv (see protein_preprocessor.py). DNA files similarly can be modified with the dna_preprocessor.py.

Dependencies: csv, argparse, shutil, pandas, numpy

Project layout (root = --folder / {project}):
  {project}/
    dna/      {case_id}.*
    rna/      (subfolders per reference manifest) / files referenced by manifest
    ch3/      (subfolders per reference manifest) / files referenced by manifest
    cnv/      (subfolders per reference manifest) / files referenced by manifest
    protein/  {case_id}.csv (or similarly named)

Manifest (case list):
  - plain text: one case ID per line (comments '#' allowed), OR
  - CSV/TSV: first column is case_id

Steps:
  --step all|dna|rna|ch3|cnv|protein
  - Running 'dna' will also chain RNA→CH3→Protein→CNV integration into DNA outputs,
    unless you pass skip flags (e.g., --skip-rna, --skip-ch3, --skip-protein, --skip-cn).

Reference:
  - By default, we look for 'reference.zip' alongside the script or in CWD.
  - You can pass --ref_zip /path/to/reference.zip OR --ref_dir /path/to/extracted/reference.
  - Inside the reference, we expect manifest files (filenames configurable via flags):
      * --rna-manifest (default: rna_manifest.tsv)
      * --ch3-manifest (default: ch3_manifest.tsv)
      * --cn-manifest  (default: cn_manifest.tsv)

Usage: 
Recommended
python multiomic_integration.py \
    --folder <project root directory> \
    --manifest <manifest file> \
    --out_dir <output directory> \
    --ref_zip <reference directory> \
    --step <all|dna|rna|ch3|cnv|protein>

Optional Single File
python multiomic_integration.py \
    --folder ./project \
    --case-id C3N-001 \
    --out_dir ./results/integration \
    --ref_dir ./reference \
    --step all

Optional
python multiomic_integration.py \
        --manifest MANIFEST_FILE \
        --input_dna_dir DNA_DIR \
        --input_rna_dir RNA_DIR \
        --rna_manifest RNA_MANIFEST \
        --input_ch3_dir CH3_DIR \
        --ch3_manifest CH3_MANIFEST \
        --input_protein_dir PROTEIN_DIR \
        --input_cn_dir CNV_DIR \
        --cn_manifest CNV_MANIFEST \
        --out_dir OUTPUT_DIR [--skip_rna] [--skip_ch3] [--skip_protein] [--skip_cn]
        

Arguments:
Required
    --folder              Project root directory
    --manifest            Case list file - formats supported: plain text (one case-ID per line), CSV/TSV (case-IDs in first column)
    --out_dir             Directory to write final integrated files
    --step                Step(s) to run [all|dna|rna|ch3|cnv|protein]

Reference Options
    --ref_dir             Directory with referrences for multiomic mapping
    --ref_zip             Path to a reference.zip archive. By default, script looks for one in current directory --ref_dir <path>

Optional
    --input_dna_dir       Directory of VEP annotated mutect CSVs named {case-ID}.csv #If not provided, --input_dna_dir defaults to dna/{case-id}.csv
    --input_rna_dir       Directory of RNA expression files
    --rna_manifest        Table linking case-IDs to RNA file paths
    --input_ch3_dir       Directory of methylation (CH3) files
    --ch3_manifest        Table linking case-IDs to CH3 file paths
    --input_protein_dir   Directory of protein annotation files
    --input_cn_dir        Directory of CNV files
    --cn_manifest         Table linking case-IDs to CNV file paths
    --skip_*              Flags to skip specific integration steps


Notes
Only the final protein files (with SNV/SNP, RNA, CH3, protein, and CNV) are written to: OUTPUT_DIR/{case-ID}.csv

"""

import os
import sys
import csv
import argparse
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import List, Optional

import pandas as pd

# ---- Utilities from pipeline ----
import logging, warnings, zipfile, shutil, numpy as np

def read_table_guess(path: Path) -> pd.DataFrame:
    """Read CSV/TSV by sniffing delimiter."""
    try:
        df = pd.read_csv(path, sep=None, engine="python")
    except Exception:
        df = pd.read_csv(path, sep="\t")
    return df

def safe_read_csv(path: str) -> Optional[pd.DataFrame]:
    try:
        return pd.read_csv(path)
    except Exception as e:
        warning(f"Could not read CSV: {path} ({e})")
        return None

def merge_with_reference(df_left: pd.DataFrame, df_right: pd.DataFrame, left_on: str, right_on: str) -> pd.DataFrame:
    return pd.merge(df_left, df_right, how='left', left_on=left_on, right_on=right_on)

def resolve_reference_dir(ref_dir_opt: Optional[str], ref_zip_opt: Optional[str]) -> Optional[Path]:
    """Resolve a directory where reference files live."""
    # Explicit directory
    if ref_dir_opt:
        ref_dir = Path(ref_dir_opt).resolve()
        if ref_dir.exists():
            info(f"Using reference directory: {ref_dir}")
            return ref_dir
        warning(f"--ref_dir not found: {ref_dir}")
    # Zip
    zcandidates = [ref_zip_opt] if ref_zip_opt else []
    zcandidates += ["reference.zip"]
    for z in zcandidates:
        if not z: 
            continue
        zp = Path(z).resolve()
        if zp.exists():
            out = zp.with_suffix("")  # ./reference
            try:
                out.mkdir(parents=True, exist_ok=True)
                with zipfile.ZipFile(zp, 'r') as zf:
                    zf.extractall(out)
                info(f"Extracted reference to: {out}")
                return out
            except Exception as e:
                warning(f"Failed to extract {zp}: {e}")
    return None


# -------------------------
# Helpers
# -------------------------

def log(msg: str) -> None:
    print(f"[multiomic_integration] {msg}")

def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)

def read_case_list(manifest: Path) -> List[str]:
    cases: List[str] = []
    with open(manifest, 'r', newline='') as f:
        sniff = f.readline()
        f.seek(0)
        if '\t' in sniff or ',' in sniff:
            # CSV/TSV: take first column
            reader = csv.reader(f, delimiter='\t' if '\t' in sniff and ',' not in sniff else ',')
            for row in reader:
                if not row: 
                    continue
                v = (row[0] or '').strip()
                if v and not v.startswith('#'):
                    cases.append(v)
        else:
            # Plain text, one per line (allow comments)
            for line in f:
                s = line.strip()
                if s and not s.startswith('#'):
                    cases.append(s.split(',')[0].strip())
    # de-duplicate, keep order
    seen = set()
    out = []
    for c in cases:
        if c not in seen:
            out.append(c); seen.add(c)
    return out

# -------------------------
# DNA loader (single-case convenience)
# -------------------------

def load_dna(folder: Path, case_id: str, input_dna_dir: Optional[str]) -> pd.DataFrame:
    # default to dna/{case-id}.csv relative to folder
    path = Path(input_dna_dir) if input_dna_dir else (folder / "dna" / f"{case_id}.csv")
    if not path.exists():
        raise FileNotFoundError(f"DNA file not found for {case_id}: {path}")
    return pd.read_csv(path)

def write_case_output(out_dir: Path, case_id: str, df: pd.DataFrame, suffix: str) -> Path:
    ensure_dir(out_dir)
    out_file = out_dir / f"{case_id}_{suffix}.csv"
    df.to_csv(out_file, index=False)
    return out_file

# -------------------------
# Placeholders for modality processors
# (Replace with your real logic if available in your project)
# -------------------------

def integrate_dna(df: pd.DataFrame) -> pd.DataFrame:
    # identity pass-through; hook for real integration
    return df

def integrate_rna(case_id: str, folder: Path, ref_dir: Optional[Path], rna_manifest_name: str, base_df: pd.DataFrame) -> Optional[pd.DataFrame]:
    if ref_dir is None:
        warning("RNA integration requested but no reference directory; skipping.")
        return base_df
    rna_manifest_path = ref_dir / rna_manifest_name
    if not rna_manifest_path.exists():
        logging.warning(f"RNA manifest not found: {rna_manifest_path}; skipping RNA.")
        return base_df
    df_rna = pd.read_table(rna_manifest_path, sep='\t', names=['Ensembl', 'Count'])
    df_rna['Ensembl'] = df_rna['Ensembl'].astype(str).str.split('.').str[0]
    merged = merge_with_reference(base_df, df_rna, 'ENSGene', 'Ensembl')
    return merged

def integrate_ch3(case_id: str, folder: Path, ref_dir: Optional[Path], ch3_manifest_name: str, base_df: pd.DataFrame) -> Optional[pd.DataFrame]:
    if ref_dir is None:
        logging.warning("CH3 integration requested but no reference directory; skipping.")
        return base_df
    ch3_manifest_path = ref_dir / ch3_manifest_name
    if not ch3_manifest_path.exists():
        logging.warning(f"CH3 manifest not found: {ch3_manifest_path}; skipping CH3.")
        return base_df
    df_ch3 = pd.read_table(ch3_manifest_path, sep='\t', names=['cg', 'beta'])
    df_ch3.rename(columns={'cg': 'IlmnID', 'beta': 'beta_val'}, inplace=True)
    merged = pd.merge(base_df, df_ch3, how='left', left_on='IlmnID', right_on='IlmnID')
    return merged

def integrate_protein(case_id: str, folder: Path, base_df: pd.DataFrame) -> Optional[pd.DataFrame]:
    p = folder / "protein" / f"{case_id}.csv"
    if p.exists():
        dfp = pd.read_csv(p)
        merged = merge_with_reference(base_df, dfp, 'To', 'NP')
        return merged
    return base_df

def integrate_cnv(case_id: str, folder: Path, ref_dir: Optional[Path], cn_manifest_name: str, base_df: pd.DataFrame) -> Optional[pd.DataFrame]:
    if ref_dir is None:
        logging.warning("CNV integration requested but no reference directory; skipping.")
        return base_df
    cn_manifest_path = ref_dir / cn_manifest_name
    if not cn_manifest_path.exists():
        logging.warning(f"CNV manifest not found: {cn_manifest_path}; skipping CNV.")
        return base_df
    df_cn = pd.read_table(cn_manifest_path), sep='\t')
    if 'copy_number' in df_cn.columns:
        df_cn['copy_number'].replace('', np.nan, inplace=True)
        df_cn.dropna(subset=['copy_number'], inplace=True)
    df_cn['EnsemblID'] = df_cn.iloc[:, 0].astype(str).str.split('.').str[0]
    cn_col = df_cn.columns[5]
    df_cn2 = df_cn[['EnsemblID', cn_col]].rename(columns={cn_col: 'copy_number'})
    merged = pd.merge(base_df, df_cn2, how='left', left_on='ENSGene', right_on='EnsemblID')
    return merged

def merge_layers(layers: List[pd.DataFrame]) -> pd.DataFrame:
    # Simplistic column-wise outer merge on shared keys if any; otherwise concat with keys
    if not layers:
        return pd.DataFrame()
    if len(layers) == 1:
        return layers[0]
    # If all have a 'gene' column, merge on that; otherwise outer-join on index
    if all('gene' in df.columns for df in layers):
        from functools import reduce
        def rename(df, tag):
            other_cols = [c for c in df.columns if c != 'gene']
            return df.rename(columns={c: f"{c}_{tag}" for c in other_cols})
        tagged = [rename(df, i) for i, df in enumerate(layers)]
        base = tagged[0]
        for nxt in tagged[1:]:
            base = base.merge(nxt, on='gene', how='outer')
        return base
    else:
        return pd.concat(layers, axis=1)

# -------------------------
# Core per-case runner
# -------------------------

def process_case(case_id: str, folder: Path, out_dir: Path, step: str,
                 input_dna_dir: Optional[str], skip_rna: bool, skip_ch3: bool, skip_protein: bool, skip_cn: bool,
                 ref_dir: Optional[Path], rna_manifest_name: str, ch3_manifest_name: str, cn_manifest_name: str) -> Path:
    log(f\"Processing {case_id} [step={step}]\" )
    outputs_dir = out_dir / case_id
    ensure_dir(outputs_dir)

    # Start with DNA
    base = pd.DataFrame()
    if step in (\"all\", \"dna\"):
        base = load_dna(folder, case_id, input_dna_dir)
        base = integrate_dna(base)

    # Chain RNA
    if step in (\"all\", \"rna\") and not skip_rna and not base.empty:
        base = integrate_rna(case_id, folder, ref_dir, rna_manifest_name, base)

    # Chain CH3
    if step in (\"all\", \"ch3\") and not skip_ch3 and not base.empty:
        base = integrate_ch3(case_id, folder, ref_dir, ch3_manifest_name, base)

    # Chain Protein
    if step in (\"all\", \"protein\") and not skip_protein and not base.empty:
        base = integrate_protein(case_id, folder, base)

    # Chain CNV
    if step in (\"all\", \"cnv\") and not skip_cn and not base.empty:
        base = integrate_cnv(case_id, folder, ref_dir, cn_manifest_name, base)

    out_file = outputs_dir / f\"{case_id}_integrated.csv\"
    base.to_csv(out_file, index=False)
    log(f\"Wrote {out_file}\")
    return out_file

# -------------------------
# CLI
# -------------------------

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Merged multiomic integration (manifest or single-case)")
    p.add_argument("--folder", required=True, help="Project root directory (contains dna/, rna/, ch3/, cnv/, protein/)")
    p.add_argument("--out_dir", required=True, help="Output directory for integrated results")
    p.add_argument("--step", default="all", choices=["all", "dna", "rna", "ch3", "cnv", "protein"], help="Which step(s) to run")

    # Multi-case via manifest OR single-case via --case-id
    p.add_argument("--manifest", help="Case list file: plain text (one per line) or CSV/TSV (case-id in first column)")
    p.add_argument("--case-id", help="Single case identifier (enables single-case mode)")

    # Single-case convenience for DNA
    p.add_argument("--input_dna_dir", default=None, help="Path to DNA CSV for single-case; defaults to dna/{case-id}.csv under --folder if not provided")

    # Reference options
    p.add_argument("--ref_dir", default=None, help="Directory with reference manifests (alternative to --ref_zip)")
    p.add_argument("--ref_zip", default=None, help="Zip containing reference manifests to extract")
    p.add_argument("--rna_manifest", default="rna_manifest.tsv", help="RNA manifest filename inside reference dir")
    p.add_argument("--ch3_manifest", default="ch3_manifest.tsv", help="CH3 manifest filename inside reference dir")
    p.add_argument("--cn_manifest", default="cn_manifest.tsv", help="CNV manifest filename inside reference dir")

    # Skip flags
    p.add_argument("--skip-rna", action="store_true", help="Skip RNA integration")
    p.add_argument("--skip-ch3", action="store_true", help="Skip CH3 integration")
    p.add_argument("--skip-protein", action="store_true", help="Skip Protein integration")
    p.add_argument("--skip-cn", action="store_true", help="Skip CNV integration")

    # Parallelism
    p.add_argument("--jobs", type=int, default=8, help="Parallel workers for multi-case mode (default: 8)")

    return p

def main():
    args = build_parser().parse_args()

    folder = Path(args.folder).resolve()
    out_dir = Path(args.out_dir).resolve()
    ensure_dir(out_dir)

    # Determine case set
    cases: List[str] = []
    if args.case_id:
        cases = [args.case_id]
    elif args.manifest:
        cases = read_case_list(Path(args.manifest).resolve())
        if not cases:
            raise SystemExit("No case IDs found in manifest")
    else:
        raise SystemExit("Provide --case-id for single-case mode or --manifest for multi-case mode.")

    # Resolve reference directory (if any)
    ref_dir = resolve_reference_dir(args.ref_dir, args.ref_zip)

    # Single-thread for single-case; parallel for multi-case
    if len(cases) == 1:
        process_case(
            cases[0], folder, out_dir, args.step,
            args.input_dna_dir,
            args.skip_rna, args.skip_ch3, args.skip_protein, args.skip_cn,
            ref_dir, args.rna_manifest, args.ch3_manifest, args.cn_manifest
        )
    else:
        jobs = max(1, int(args.jobs))
        with ProcessPoolExecutor(max_workers=jobs) as ex:
            futs = {
                ex.submit(
                    process_case, cid, folder, out_dir, args.step,
                    args.input_dna_dir,
                    args.skip_rna, args.skip_ch3, args.skip_protein, args.skip_cn,
                    ref_dir, args.rna_manifest, args.ch3_manifest, args.cn_manifest
                ): cid for cid in cases
            }
            for fut in as_completed(futs):
                cid = futs[fut]
                try:
                    fut.result()
                except Exception as e:
                    log(f"ERROR [{cid}]: {e}")

if __name__ == "__main__":
    main()
