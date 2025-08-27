#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
DNA PREPROCESSOR

This script creates per-case DNA CSVs under 'dna/{Case-ID}.csv' that can be used in multiomic integration. It reads a GDC-like 'dna_manifest.tsv' and parses VCF/VCF.GZ files for each case, using a VCF parser to select for variants that pass standard filters (PASS and alt_allele_in_normal). 

Dependencies: pandas

Usage:
python dna_preprocessor.py \
  --in_dir <input directory> \
  --manifest <gdc-manifest tsv> \
  --out_dir <output directory>

(Optional)\
  --make-simplified \
  --simplified <case_ids.txt> \

Arguments:
Required
   --manifest           GDC-like TSV/CSV with at least Case ID, File Name (File ID optional)
   --in_dir             Input directory that contains raw data, format <in_dir>/<File ID>/<File Name>
   --out_dir            Output directory that will contain <dna/<Case-ID>.csv that can be used in multiomic integration

General
   --max-records N      Cap parsed VCF rows per case (for smoke tests)
   --jobs               Controls parallel execution, if not provided, script uses min(8, CPU count)

Make/list Case-IDS:
   --make-simplified    Provides unique Case-IDs derived from --manifest
   --simplified         Path to write the Case-ID list (default: <out_dir>/case_ids.txt)

Notes
1. Case-ID normalization: uses first token before a comma (e.g., "case-01, C3N-04155" -> case-01)
2. VCF parser: minimal; captures core fields; genotype fields are not parsed in the main CSVs.
3. Duplicates: first row per Case-ID in the manifest "wins".

Troubleshooting
1. “DNA source file not found …” → Check that File Name and (if used) File ID match your folder/dna or folder/vep paths, or provide an absolute path in File Name.
2. No CSVs produced → Verify Case ID and File Name columns exist in your manifest and that files are .vcf or .vcf.gz.

"""
from concurrent.futures import ProcessPoolExecutor, as_completed
import os
import argparse
import logging
import gzip
import zipfile
import shutil
import tempfile
from contextlib import contextmanager, nullcontext
from typing import Optional, List, Tuple
from pathlib import Path
import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")


# --------------------------------------------------
# Core VCF parsing
# --------------------------------------------------

def parse_vcf_to_df(vcf_path: Path, max_records: Optional[int] = None) -> pd.DataFrame:
    """
    Minimal VCF parser -> DataFrame with columns: CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO.
    """
    rows: List[dict] = []
    opener = gzip.open if str(vcf_path).endswith(".gz") else open
    with opener(vcf_path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            chrom, pos, vid, ref, alt, qual, flt, info = parts[:8]
            rows.append({
                "CHROM": chrom,
                "POS": int(pos) if pos.isdigit() else pos,
                "ID": vid,
                "REF": ref,
                "ALT": alt,
                "QUAL": qual,
                "FILTER": flt,
                "INFO": info,
            })
            if max_records is not None and len(rows) >= max_records:
                break
    return pd.DataFrame(rows, columns=["CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO"])


# --------------------------------------------------
# File helpers
# --------------------------------------------------

def _norm_case_id(val: str) -> str:
    s = str(val).strip().strip('"').strip("'")
    if "," in s:
        s = s.split(",", 1)[0].strip()
    return s

def emit_simplified_case_list(manifest_path: Path, out_path: Path) -> Path:
    """
    Extract unique Case-IDs from a dna_manifest.tsv (or similar) and write one per line.
    - Reads CSV/TSV by sniffing delimiter (fallback TSV)
    - Prefers named columns: Case ID / Case-ID / CaseID / ID
    - Else guesses column index (6th if present, else 1st)
    - Normalizes (strip quotes, token before comma), de-dupes preserving order
    """
    manifest_path = Path(manifest_path)
    out_path = Path(out_path)

    try:
        df = pd.read_csv(manifest_path, sep=None, engine="python", dtype=str)
    except Exception:
        df = pd.read_csv(manifest_path, sep="\t", dtype=str)
    df = df.fillna("")

    cols = {c.lower().strip(): c for c in df.columns}
    case_col: Optional[str] = None
    for key in ("case id", "case-id", "caseid", "id"):
        if key in cols:
            case_col = cols[key]
            break

    if case_col is None:
        idx = 5 if df.shape[1] > 5 else 0
        series = df.iloc[:, idx].astype(str)
    else:
        series = df[case_col].astype(str)

    ids_norm: List[str] = []
    seen = set()
    for x in series.tolist():
        x = x.strip()
        if not x:
            continue
        nx = _norm_case_id(x)
        if nx and nx not in seen:
            seen.add(nx)
            ids_norm.append(nx)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(ids_norm) + ("\n" if ids_norm else ""))
    return out_path

def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)

def read_table_guess(path: Path) -> pd.DataFrame:
    """Read CSV/TSV by sniffing delimiter (fallback to TSV)."""
    try:
        return pd.read_csv(path, sep=None, engine="python")
    except Exception:
        return pd.read_csv(path, sep="\t")

def resolve_dna_file(project_root: Path, file_id: str, file_name: str) -> Path:
    """
    Try common layouts:
      - project_root/dna/{File ID}/{File Name}
      - project_root/vep/{File ID}/{File Name}
      - direct path
      - search by basename anywhere under project_root
      - handle .vcf.gz vs .vcf swaps
      - fallback: search by partial File ID
    """
    root = Path(project_root)

    candidate = root / file_name
    if candidate.exists():
        return candidate

    base = Path(file_name).name
    hit = next(root.rglob(base), None)
    if hit is not None:
        return hit

    if base.endswith(".vcf.gz"):
        alt = base[:-3]
        hit = next(root.rglob(alt), None)
        if hit is not None:
            return hit
    if base.endswith(".vcf"):
        alt = base + ".gz"
        hit = next(root.rglob(alt), None)
        if hit is not None:
            return hit

    fid = (file_id or "").strip()
    if fid:
        hit = next((p for p in root.rglob("*") if p.is_file() and fid in p.name), None)
        if hit is not None:
            return hit

    raise FileNotFoundError(f"Could not resolve DNA file for '{file_name}' (id='{file_id}') under {root}")

# ---- compression helpers ----

COMPRESSED_EXTS = (".gz", ".zip")
TEXT_LIKE = (".vcf", ".tsv", ".txt", ".maf", ".csv")

def is_zip(p: Path) -> bool: return p.suffix.lower() == ".zip"
def is_compressed(p: Path) -> bool: return p.suffix.lower() in COMPRESSED_EXTS

def _strip_gz_suffix(p: Path) -> Path:
    return p.with_suffix("") if p.suffix.lower() == ".gz" else p

@contextmanager
def _scratch(base_dir: Optional[str]) -> Path:
    """
    Yields a scratch directory Path. If base_dir is provided, re-use it; otherwise
    create a temp directory and clean it up automatically.
    """
    if base_dir:
        d = Path(base_dir).resolve()
        d.mkdir(parents=True, exist_ok=True)
        yield d
    else:
        with tempfile.TemporaryDirectory(prefix="dna_preproc_") as td:
            yield Path(td)

def materialize_input(p: Path, scratch_dir: Path) -> Path:
    """
    Ensure an input path is a real, readable file on disk:
      - If uncompressed: return as-is.
      - If .gz: decompress to scratch and return decompressed path.
      - If .zip: extract to scratch and return the most likely text-like file.
    """
    if not p.exists():
        raise FileNotFoundError(f"Input not found: {p}")

    if not is_compressed(p):
        return p

    if is_zip(p):
        extract_root = scratch_dir / (p.stem + "_unzipped")
        extract_root.mkdir(parents=True, exist_ok=True)
        with zipfile.ZipFile(p, "r") as zf:
            zf.extractall(extract_root)
        candidates: List[Path] = []
        for ext in TEXT_LIKE:
            candidates.extend(extract_root.rglob(f"*{ext}"))
        if candidates:
            candidates.sort(key=lambda x: len(str(x)))
            return candidates[0]
        all_files = [q for q in extract_root.rglob("*") if q.is_file()]
        if all_files:
            all_files.sort(key=lambda x: len(str(x)))
            return all_files[0]
        raise FileNotFoundError(f"Zip archive had no files: {p}")

    # .gz
    out_path = scratch_dir / _strip_gz_suffix(p).name
    with gzip.open(p, "rb") as src, open(out_path, "wb") as dst:
        shutil.copyfileobj(src, dst)
    return out_path


# --------------------------------------------------
# Parallel worker
# --------------------------------------------------

def _process_one_case(
    case_id: str,
    file_id: str,
    file_name: str,
    project_root: str,
    out_root: str,
    max_records: Optional[int],
    unzip_inputs: bool,
    base_scratch_dir: Optional[str],
) -> Tuple[str, str, int, int, str]:
    try:
        safe_case = case_id.replace("/", "_")
        project_root_p = Path(project_root)
        out_root_p = Path(out_root)
        dna_dir = out_root_p / "dna"
        ensure_dir(dna_dir)

        src = resolve_dna_file(project_root_p, file_id, file_name)

        # Materialize into a per-case scratch subfolder when requested
        eff_src = src
        if unzip_inputs:
            case_scratch = (Path(base_scratch_dir) / safe_case) if base_scratch_dir else (out_root_p / f"_scratch_{os.getpid()}_{safe_case}")
            case_scratch.mkdir(parents=True, exist_ok=True)
            eff_src = materialize_input(src, case_scratch)

        vdf = parse_vcf_to_df(eff_src, max_records=max_records)

        # Add Case-ID, filter PASS ∪ alt_allele_in_normal
        vdf.insert(0, "Case-ID", case_id)
        fvals = vdf["FILTER"].astype(str).str.strip()
        is_pass = fvals.str.upper().eq("PASS")
        is_alt  = fvals.str.contains("alt_allele_in_normal", case=False, na=False)
        vdf_keep = vdf[is_pass | is_alt].copy()

        if vdf_keep.empty:
            return ("skip", str(dna_dir / f"{safe_case}.csv"), 0, len(vdf), "No PASS/alt_allele_in_normal rows")

        out_path = dna_dir / f"{safe_case}.csv"
        vdf_keep.to_csv(out_path, index=False)
        return ("ok", str(out_path), len(vdf_keep), len(vdf), "")

    except Exception as e:
        return ("err", "", 0, 0, f"{type(e).__name__}: {e}")


# --------------------------------------------------
# CLI
# --------------------------------------------------

def parse_args():
    default_workers = min(8, os.cpu_count() or 1)
    p = argparse.ArgumentParser(description="Create per-case DNA CSVs (dna/{Case-ID}.csv) from VCFs.")
    p.add_argument("--in_dir", required=True, help="Project root containing dna/ or vep/ subfolders")
    p.add_argument("--manifest", required=True, help="Path to dna_manifest.tsv (GDC-like)")
    p.add_argument("--out_dir", required=True, help="Output directory root (CSV written under out_dir/dna/)")
    p.add_argument("--max-records", type=int, default=None, help="Optional cap on parsed VCF rows per case")
    p.add_argument("--make-simplified", action="store_true", help="Emit a Case-ID list derived from --manifest")
    p.add_argument("--simplified", help="Path to write the Case-ID list (default: <out_dir>/case_ids.txt)")
    p.add_argument("--jobs", type=int, default=default_workers, help=f"Parallel workers (default={default_workers})")
    p.add_argument("--unzip-inputs", action="store_true", help="Materialize .gz/.zip inputs into scratch before parsing")
    p.add_argument("--scratch-dir", dest="scratch_dir", default=None, help="Base scratch directory (per-case subdirs are created within)")
    return p.parse_args()


# --------------------------------------------------
# Main
# --------------------------------------------------

def main():
    args = parse_args()
    project_root = Path(args.in_dir).resolve()
    out_root = Path(args.out_dir).resolve()
    ensure_dir(out_root / "dna")

    # Read manifest
    df = read_table_guess(Path(args.manifest))

    # Optional: emit simplified list (robust)
    simp_out = Path(args.simplified) if args.simplified else (out_root / "case_ids.txt")
    if args.make_simplified:
        outp = emit_simplified_case_list(Path(args.manifest), simp_out)
        logging.info(f"Wrote simplified Case-ID list: {outp}")

    # column detection
    cols = {c.lower().strip(): c for c in df.columns}
    case_col = next((cols[k] for k in ["case id","case-id","caseid"] if k in cols), None)
    if case_col is None:
        raise SystemExit("dna_manifest.tsv must include a 'Case ID' column")
    file_id_col = next((cols[k] for k in ["file id","file-id","fileid"] if k in cols), None)
    file_name_col = next((cols[k] for k in ["file name","file-name","filename"] if k in cols), None)
    if file_name_col is None:
        raise SystemExit("dna_manifest.tsv must include a 'File Name' column")

    # Build task list (first occurrence per Case-ID wins)
    tasks: List[Tuple[str, str, str]] = []
    seen = set()
    for _, row in df.iterrows():
        case_id = str(row[case_col]).split(",")[0].strip().strip('"')
        if not case_id or case_id in seen:
            continue
        seen.add(case_id)
        file_id = str(row[file_id_col]).strip() if file_id_col and pd.notna(row[file_id_col]) else ""
        file_name = str(row[file_name_col]).strip()
        if not file_name:
            logging.warning(f"[{case_id}] Missing 'File Name'; skipping")
            continue
        tasks.append((case_id, file_id, file_name))

    if not tasks:
        logging.warning("No cases queued from manifest.")
        return

    # Prepare a base scratch directory (shared), then workers use per-case subdirs
    scratch_cm = _scratch(args.scratch_dir) if args.unzip_inputs else nullcontext(Path())
    with scratch_cm as scratch_base:
        base_scratch_str = str(scratch_base) if args.unzip_inputs else None

        ok = skip = err = 0
        with ProcessPoolExecutor(max_workers=args.jobs) as ex:
            futures = [
                ex.submit(
                    _process_one_case,
                    cid, fid, fname,
                    str(project_root), str(out_root),
                    args.max_records,
                    args.unzip_inputs,
                    base_scratch_str,
                )
                for cid, fid, fname in tasks
            ]
            for fut in as_completed(futures):
                status, outp, kept, total, msg = fut.result()
                if status == "ok":
                    ok += 1
                    logging.info(f"[ok] wrote {outp} ({kept}/{total})")
                elif status == "skip":
                    skip += 1
                    logging.warning(f"[skip] {msg}")
                else:
                    err += 1
                    logging.error(f"[err] {msg}")

    logging.info(f"Done. wrote={ok}, skipped={skip}, errors={err}. Output: {out_root/'dna'}")


if __name__ == "__main__":
    main()
