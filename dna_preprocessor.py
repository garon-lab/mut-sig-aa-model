#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
DNA PREPROCESSOR

This script creates per-case DNA CSVs under 'dna/{Case-ID}.csv' that can be used in multiomic integration. It reads a GDC-like 'dna_manifest.tsv' and parses VCF/VCF.GZ files for each case. Our manuscript selects for single nucleotide variants (SNVs) and sningle nucleotide polymorphisms (SNPs), which is available pre-built in optional steps, but can be adapted as needed by the user. 

This script processes VEP annotated mutect files by:
1. Using a VCF parser to select for variants of interest (e.g., SNV/SNP).
2. Summarizing SNP/SNV counts (optional).
3. Extracting DNA mutational signatures (optional).
4. Extracting amino acid substitutions from SNV/SNPs (optional).
5. Displaying in amino acid substitutions in matrix format for further comparison/analysis (optional).

Dependencies: pandas

Usage:
python dna_preprocessor.py \
  --folder <input directory> \
  --manifest <gdc-manifest tsv> \
  --out_dir <output directory>

(Optional)\
  --make-simplified \
  --preprocess-mutect \
  --vcf-folder <vcf directory> \
  --simplified <case_ids.txt> \
  --write-signatures --signature-label <type> defaults dna> \
  --extract-mutations <type, e.g., snp|snv> --write-matrices 

Arguments:
(Required)\
   --manifest           GDC-like TSV/CSV with at least Case ID, File Name (File ID optional)\
   --folder             Input directory that contains raw data, format <folder>/dna/<File ID>/<File Name>\
   --out_dir            Output directory that will contain <dna/<Case-ID>.csv that can be used in multiomic integration\

General:\
   --max-records N      Cap parsed VCF rows per case (for smoke tests)\
   --jobs               Controls parallel execution, if not provided, script uses min(8, CPU count)

Make/list Case-IDS:\
   --make-simplified    Provides unique Case-IDs derived from --manifest\
   --simplified         Path to write the Case-ID list (default: <out_dir>/case_ids.txt)\
   
Preprocess: \
   --preprocess-mutect  Flag for extended analysis, strips '##' headers and writes prep/<Case-ID>.txt\
   --vcf-folder         Where per-case VCFs live (default <folder/dna>)

Analytics (require --simplified file listing Case-IDs):\
   --simplified FILE              Path to case_ids.txt if it has been previously made, should have one Case-ID per line (no header)\
   --summarize-variants           Write SNP/SNV counts to <out.dir>/summary.csv\
   --write-signatures             Write <out_dir>/<label>-signature.csv\
   --signature-label L            Label for signature file prefix (default: dna)\
   --extract-mutations {snp|snv}  Extracts ST/END AA pairs to <out_dir>/<type>/<Case-ID>.csv\
   --write-matrices               Writes 21 x 21 amino acid matrices to <out_dir>/<type>/matrices/<Case-IDs>.csv\

Notes
1. Case-ID normalization: uses first token before a comma (e.g., "case-01, C3N-04155" -> case-01)
2. VCF parser: minimal; catpures core fields; genotype fileds are not parsed in the main CSVs.
3. Analytics: expects MuTect-style flast TSVs from --preprocess-mutect with columns mapped to ['#CHROM', 'POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','NORMAL','TUMOR']
4. Filters: SNP - FILTER contains alt and INFO contains missense, SNV - FILTER contains PASS and INFO contains missense
5. Amino acid (AA) parsing: extract-mutations uses the 16th INFO pipe-filed (info.split('|')[15]) and takes its first/last char as AA start/end. Adjust if your annotation format differs.
6. Duplicates: first row per Case-ID in the manifest "wins".

Troubleshooting
1. “DNA source file not found …” → Check that File Name and (if used) File ID match your folder/dna or folder/vep paths, or provide an absolute path in File Name.
2. No CSVs produced → Verify Case ID and File Name columns exist in your manifest and that files are .vcf or .vcf.gz.
3. Optional steps say file not found → Run --preprocess-mutect first to populate <out_dir>/prep/*.txt or ensure those files already exist.
4. AA field index errors → Your annotation pipeline may place AA changes in a different INFO index. Update the parsing logic accordingly.

"""

import argparse
import logging
import gzip
import zipfile
import shutil
import tempfile
from contextlib import contextmanager
from contextlib import nullcontext
from typing import Optional, List, Tuple
from pathlib import Path
import pandas as pd
import gzip

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

# Amino-acid list for matrices (single-letter codes)
AA_LIST = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']

def summarize_variants(simplified: Path, out_dir: Path):
    """
    Optional step: summarize SNP/SNV counts per Case-ID from preprocessed VCFs under out_dir/prep/*.txt
    Writes: out_dir/summary.csv with header: Case-ID,SNP,SNV
    """
    import pandas as _pd
    simplified = Path(simplified)
    out_dir = Path(out_dir)
    df = _pd.read_csv(simplified, header=None)
    ids = df.iloc[:, 0]
    summary_file = out_dir / "summary.csv"
    summary_file.parent.mkdir(parents=True, exist_ok=True)
    with summary_file.open('w') as f:
        f.write("Case-ID,SNP,SNV\n")

    for case_id in ids:
        try:
            filepath = out_dir / "prep" / f"{case_id}.txt"
            vcf_df = _pd.read_table(filepath, sep='\t', header=None)
            vcf_df.columns = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'NORMAL', 'TUMOR']
            snp_count = vcf_df[(vcf_df['FILTER'].str.contains("alt", na=False)) & (vcf_df['INFO'].str.contains("missense", na=False))].shape[0]
            snv_count = vcf_df[(vcf_df['FILTER'].str.contains("PASS", na=False)) & (vcf_df['INFO'].str.contains("missense", na=False))].shape[0]
            with summary_file.open('a') as f:
                f.write(f"{case_id},{snp_count},{snv_count}\n")
        except Exception as e:
            logging.warning(f"Skipping {case_id}: {e}")

def write_signatures(prep_dir: Path, simplified: Path, out_dir: Path, label: str):
    """
    Optional step: compute simple mutation 'signatures' per case from preprocessed VCFs.
    Writes: out_dir/{label}-signature.csv with header: Case-ID,SUM,CTGA,CAGT,GCCG,ATTA,AGTC,ACTG
    """
    import pandas as _pd
    prep_dir = Path(prep_dir)
    simplified = Path(simplified)
    out_dir = Path(out_dir)
    sample_ids = _pd.read_csv(simplified, header=None).iloc[:, 0]
    out_file = out_dir / f"{label}-signature.csv"
    header = ['Case-ID', 'SUM', 'CTGA', 'CAGT', 'GCCG', 'ATTA', 'AGTC', 'ACTG']
    all_rows = []
    for case_id in sample_ids:
        file_path = prep_dir / f"{case_id}.txt"
        if not file_path.exists():
            logging.warning(f"{case_id}: file not found at {file_path}")
            continue
        try:
            df = _pd.read_table(file_path)
            df['sA'] = df.iloc[:, 3].str.contains('A', na=False)
            df['sC'] = df.iloc[:, 3].str.contains('C', na=False)
            df['sG'] = df.iloc[:, 3].str.contains('G', na=False)
            df['sT'] = df.iloc[:, 3].str.contains('T', na=False)
            df['eA'] = df.iloc[:, 4].str.contains('A', na=False)
            df['eC'] = df.iloc[:, 4].str.contains('C', na=False)
            df['eG'] = df.iloc[:, 4].str.contains('G', na=False)
            df['eT'] = df.iloc[:, 4].str.contains('T', na=False)
            df['NM'] = df.iloc[:, 6].str.strip().str[:3]
            df1 = df[df['NM'].str.contains('alt', na=False)]
            mutations = {
                'AC': df1['sA'] & df1['eC'],
                'AG': df1['sA'] & df1['eG'],
                'AT': df1['sA'] & df1['eT'],
                'CA': df1['sC'] & df1['eA'],
                'CG': df1['sC'] & df1['eG'],
                'CT': df1['sC'] & df1['eT'],
                'GA': df1['sG'] & df1['eA'],
                'GC': df1['sG'] & df1['eC'],
                'GT': df1['sG'] & df1['eT'],
                'TA': df1['sT'] & df1['eA'],
                'TC': df1['sT'] & df1['eC'],
                'TG': df1['sT'] & df1['eG'],
            }
            counts = {
                'CTGA': (mutations['CT'] | mutations['GA']).sum(),
                'CAGT': (mutations['CA'] | mutations['GT']).sum(),
                'GCCG': (mutations['CG'] | mutations['GC']).sum(),
                'ATTA': (mutations['AT'] | mutations['TA']).sum(),
                'AGTC': (mutations['AG'] | mutations['TC']).sum(),
                'ACTG': (mutations['AC'] | mutations['TG']).sum(),
            }
            total = sum(counts.values())
            row = [case_id, total] + [counts[key] for key in ['CTGA', 'CAGT', 'GCCG', 'ATTA', 'AGTC', 'ACTG']]
            all_rows.append(row)
        except Exception as e:
            logging.warning(f"Failed to process {case_id}: {e}")
    df_out = _pd.DataFrame(all_rows, columns=header)
    df_out.to_csv(out_file, index=False)

def extract_mutations(prep_dir: Path, out_dir: Path, simplified: Path, mutation_type: str):
    """
    Optional step: extract amino-acid pairs per mutation type (snp|snv) from preprocessed VCFs.
    Writes: out_dir/{mtype}/{Case-ID}-{mtype}.csv with columns: ST, END, #CHROM, TUMOR
    """
    assert mutation_type in ["snp", "snv"]
    import pandas as _pd
    prep_dir = Path(prep_dir)
    out_dir = Path(out_dir)
    simplified = Path(simplified)
    df = _pd.read_csv(simplified, header=None)
    ids = df.iloc[:, 0]
    (out_dir / mutation_type).mkdir(parents=True, exist_ok=True)
    for case_id in ids:
        try:
            infile = prep_dir / f"{case_id}.txt"
            vcf_df = _pd.read_table(infile, sep='\t', header=None)
            vcf_df.columns = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'NORMAL', 'TUMOR']
            if mutation_type == "snp":
                filtered = vcf_df[vcf_df['FILTER'].str.contains("alt", na=False)]
            else:
                filtered = vcf_df[vcf_df['FILTER'].str.contains("PASS", na=False)]
            info = filtered['INFO'].str.split('|', expand=True)
            if 15 not in info.columns:
                logging.warning(f"Case {case_id}: INFO field has fewer than 16 fields - skipping")
                continue
            aa_df = _pd.DataFrame({
                'ST': info[15].str[0],
                'END': info[15].str[-1],
                '#CHROM': filtered['#CHROM'],
                'TUMOR': filtered['TUMOR']
            })
            outfile = out_dir / mutation_type / f"{case_id}-{mutation_type}.csv"
            aa_df.to_csv(outfile, index=False)
        except Exception as e:
            logging.warning(f"Failed to extract {mutation_type} for {case_id}: {e}")

def generate_aa_matrix(df):
    matrix = df.groupby(['ST', 'END']).size().unstack(fill_value=0)
    matrix = matrix.reindex(index=AA_LIST + ['*'], columns=AA_LIST + ['*'], fill_value=0)
    return matrix

def write_matrices(out_dir: Path, simplified: Path):
    """
    Optional step: write per-case amino-acid matrices under out_dir/snp/matrices and out_dir/snv/matrices.
    """
    import pandas as _pd
    out_dir = Path(out_dir)
    simplified = Path(simplified)
    df = _pd.read_csv(simplified, header=None)
    ids = df.iloc[:, 0]
    for case_id in ids:
        for mtype in ['snp', 'snv']:
            try:
                file = out_dir / mtype / f"{case_id}-{mtype}.csv"
                if not file.exists():
                    continue
                df_aa = _pd.read_csv(file)
                matrix = generate_aa_matrix(df_aa)
                out_path = out_dir / mtype / "matrices" / f"{case_id}.csv"
                out_path.parent.mkdir(parents=True, exist_ok=True)
                matrix.to_csv(out_path)
            except Exception as e:
                logging.warning(f"Matrix generation failed for {case_id}-{mtype}: {e}")
              
def _find_vcf_for_case(folder: Path, case_id: str) -> Optional[Path]:
    """
    Try common layouts:
      - <folder>/{case_id}.vcf
      - <folder>/{case_id}.vcf.gz
      - any depth: **/*{case_id}*.vcf or **/*{case_id}*.vcf.gz
    Returns a Path or None if not found.
    """
    folder = Path(folder)
    # Fast checks first
    direct = [folder / f"{case_id}.vcf", folder / f"{case_id}.vcf.gz"]
    for p in direct:
        if p.exists():
            return p
    # Recursive fallback: relax to contains case_id
    # (adjust if you need stricter matching)
    for pat in (f"*{case_id}*.vcf", f"*{case_id}*.vcf.gz"):
        hit = next(folder.rglob(pat), None)
        if hit is not None:
            return hit
    return None

def preprocess_mutect(folder: Path, manifest_file: Path, out_dir: Path):
    """
    Optional step: strip '##' lines from per-case VCFs named {Case-ID}.vcf[.gz] to prep/{Case-ID}.txt
    - folder: directory containing {Case-ID}.vcf or {Case-ID}.vcf.gz
    - manifest_file: dna_manifest.tsv (GDC-like) or a simple list; we will extract Case IDs
    - out_dir: root output directory (prep files written to out_dir/prep)
    """
    import pandas as _pd, gzip as _gzip
    folder = Path(folder)
    out_dir = Path(out_dir)
    prep_dir = out_dir / "prep"
    prep_dir.mkdir(parents=True, exist_ok=True)
    # Read manifest (dna_manifest.tsv) and get a first-column list if no header
    try:
        df = _pd.read_csv(manifest_file, sep=None, engine="python")
        cols = {c.lower().strip(): c for c in df.columns}
        if "case id" in cols: 
            ids = df[cols["case id"]].astype(str)
        elif "case-id" in cols:
            ids = df[cols["case-id"]].astype(str)
        elif "caseid" in cols:
            ids = df[cols["caseid"]].astype(str)
        else:
            # headerless fallback: take column 0 or 5 if looks like GDC export
            pick_col = 0 if 0 < df.shape[1] else None
            if df.shape[1] > 5:
                pick_col = 5
            ids = df.iloc[:, pick_col].astype(str)
    except Exception:
        # Plain text list
        ids = [line.strip() for line in Path(manifest_file).read_text().splitlines() if line.strip()]

    for case_id in ids:
      case_id = str(case_id).split(",")[0].strip().strip('"')
      vcf_path = _find_vcf_for_case(folder, case_id)
      target = (Path(out_dir) / "prep" / f"{case_id}.txt")
      try:
          if vcf_path is None:
              logging.warning(f"{case_id}: VCF not found under {folder}")
              continue
          if vcf_path.suffix == ".gz":
              with _gzip.open(vcf_path, "rt", encoding="utf-8", errors="replace") as fin, open(target, "w") as fout:
                  for line in fin:
                      if not line.startswith("##"):
                          fout.write(line)
          else:
              with open(vcf_path, "rt", encoding="utf-8", errors="replace") as fin, open(target, "w") as fout:
                  for line in fin:
                      if not line.startswith("##"):
                          fout.write(line)
      except Exception as e:
          logging.warning(f"{case_id}: failed to preprocess VCF: {e}")


# ---- Helpers ----


def _norm_case_id(val: str) -> str:
    s = str(val).strip().strip('"').strip("'")
    if "," in s:
        s = s.split(",")[0].strip()
    return s

def emit_simplified_case_list(manifest_path: Path, out_path: Path) -> Path:
    """
    Extract unique Case-IDs from a dna_manifest.tsv (or similar) and write one per line.
    If the file appears headerless, we fall back to a reasonable column guess.
    """
    import pandas as _pd
    manifest_path = Path(manifest_path)
    out_path = Path(out_path)
    try:
        df = _pd.read_csv(manifest_path, sep=None, engine="python")
    except Exception:
        df = _pd.read_csv(manifest_path, sep="\t")

    # Try to find a Case ID column (case-insensitive)
    cols = {c.lower().strip(): c for c in df.columns}
    case_col = None
    for key in ["case id", "case-id", "caseid", "id"]:
        if key in cols:
            case_col = cols[key]
            break

    if case_col is None:
        # Headerless or unexpected headers -> guess a column index
        # Prefer column 5 (Case ID in many GDC exports), else column 0
        try:
            if df.shape[1] > 5:
                ids = df.iloc[:, 5].astype(str)
            else:
                ids = df.iloc[:, 0].astype(str)
        except Exception:
            # Last resort: treat file as plain text list
            ids = _pd.Series([line.strip() for line in manifest_path.read_text().splitlines() if line.strip()])
    else:
        ids = df[case_col].astype(str)

    # Normalize and uniquify
    ids = [_norm_case_id(x) for x in ids if str(x).strip()]
    uniq = sorted(dict.fromkeys(ids))

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(uniq) + "\n")
    return out_path

COMPRESSED_EXTS = (".gz", ".zip")
TEXT_LIKE = (".vcf", ".tsv", ".txt", ".maf", ".csv")

def _strip_gz_suffix(p: Path) -> Path:
    # Remove a single .gz suffix
    if p.suffix == ".gz":
        return p.with_suffix("")  # e.g., foo.vcf.gz -> foo.vcf
    return p

def is_compressed(p: Path) -> bool:
    s = p.suffix.lower()
    return s in COMPRESSED_EXTS

def is_zip(p: Path) -> bool:
    return p.suffix.lower() == ".zip"

@contextmanager
def _scratch(args) -> Path:
    """
    Yields a scratch directory Path. If --scratch-dir is provided, re-use it.
    Otherwise, create a TemporaryDirectory and clean it up automatically.
    """
    if args.scratch_dir:
        d = Path(args.scratch_dir).resolve()
        d.mkdir(parents=True, exist_ok=True)
        yield d
    else:
        with tempfile.TemporaryDirectory(prefix="dna_preproc_") as td:
            yield Path(td)

def materialize_input(p: Path, scratch_dir: Path) -> Path:
    """
    Ensure an input path is a real, readable file on disk.
    - If uncompressed: return as-is (after existence check).
    - If .gz: decompress to scratch and return the decompressed path.
    - If .zip: extract to a subdir and return the most likely file (heuristic).
    """
    if not p.exists():
        # Some manifests store relative paths — try to be forgiving
        # (Callers should also join a project root if needed)
        raise FileNotFoundError(f"Input not found: {p}")

    if not is_compressed(p):
        return p

    if is_zip(p):
        extract_root = scratch_dir / (p.stem + "_unzipped")
        extract_root.mkdir(parents=True, exist_ok=True)
        with zipfile.ZipFile(p, "r") as zf:
            zf.extractall(extract_root)
        # Heuristic: pick a likely text-like file if present
        candidates: List[Path] = []
        for ext in TEXT_LIKE:
            candidates.extend(extract_root.rglob(f"*{ext}"))
        if candidates:
            # Prefer shortest path / “closest match”
            candidates.sort(key=lambda x: len(str(x)))
            return candidates[0]
        # Fallback: first file in tree
        all_files = [q for q in extract_root.rglob("*") if q.is_file()]
        if all_files:
            all_files.sort(key=lambda x: len(str(x)))
            return all_files[0]
        raise FileNotFoundError(f"Zip archive had no files: {p}")

    # .gz case: stream-decompress to scratch
    out_path = scratch_dir / _strip_gz_suffix(p).name
    with gzip.open(p, "rb") as src, open(out_path, "wb") as dst:
        shutil.copyfileobj(src, dst)
    return out_path

def read_table_guess(path: Path) -> pd.DataFrame:
    """Read CSV/TSV by sniffing delimiter (fallback to TSV)."""
    try:
        return pd.read_csv(path, sep=None, engine="python")
    except Exception:
        return pd.read_csv(path, sep="\t")

def parse_vcf_to_df(vcf_path: Path, max_records: int = None) -> pd.DataFrame:
    """
    Minimal VCF parser -> DataFrame.
    Captures CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO (raw).

    If FORMAT/sample columns exist, they are not parsed here (kept simple).
    """
    rows = []
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
                "INFO": info
            })
            if max_records is not None and len(rows) >= max_records:
                break
    return pd.DataFrame(rows, columns=["CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO"])

def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)

def resolve_dna_file(project_root: Path, file_id: str, file_name: str) -> Path:
    """
    Try several common layouts:
      1) {project_root}/dna/{File ID}/{File Name}
      2) {project_root}/vep/{File ID}/{File Name}
      3) Absolute path in file_name if it exists
    """
    root = Path(project_root)

    # 1) direct exact path if file_name contains dirs
    candidate = root / file_name
    if candidate.exists():
        return candidate

    # 2) search by basename anywhere under root
    base = Path(file_name).name
    hit = next(root.rglob(base), None)
    if hit is not None:
        return hit

    # 3) if it's .vcf.gz in the manifest, try plain .vcf in case it was decompressed
    if base.endswith(".vcf.gz"):
        alt = base[:-3]  # drop .gz -> .vcf
        hit = next(root.rglob(alt), None)
        if hit is not None:
            return hit

    # 4) if it's .vcf in the manifest, try .vcf.gz
    if base.endswith(".vcf"):
        alt = base + ".gz"
        hit = next(root.rglob(alt), None)
        if hit is not None:
            return hit

    # 5) last-ditch: search by (partial) file_id token if present
    fid = (file_id or "").strip()
    if fid:
        hit = next((p for p in root.rglob("*") if p.is_file() and fid in p.name), None)
        if hit is not None:
            return hit

    raise FileNotFoundError(f"Could not resolve DNA file for '{file_name}' (id='{file_id}') under {root}")

# ---- Core ----

def parse_args():
    p = argparse.ArgumentParser(description="Create per-case DNA CSVs (dna/{Case-ID}.csv) from a dna_manifest.tsv using a VCF parser.")
    p.add_argument("--folder", required=True, help="Project root containing dna/ or vep/ subfolders")
    p.add_argument("--manifest", required=True, help="Path to dna_manifest.tsv (GDC-like)")
    p.add_argument("--out_dir", required=True, help="Output directory root (CSV written under out_dir/dna/)")
    p.add_argument("--unzip-inputs", action="store_true", help="If set, automatically materialize compressed inputs (.gz, .zip) into a scratch directory before processing.")
    p.add_argument("--scratch-dir",default=None, help="Optional directory to use for temporary extracted files. Defaults to a new Temporary Directory each run.")
    p.add_argument("--max-records", type=int, default=None, help="Optional cap on parsed VCF records per case (for testing)")
    p.add_argument("--make-simplified", action="store_true", help="Emit a Case-ID list derived from --manifest")
    p.add_argument("--simplified", help="Path to write the Case-ID list (default: <out_dir>/case_ids.txt)")
    p.add_argument("--preprocess-mutect", dest="preprocess_mutect", action="store_true", help="Preprocess Mutect VCFs prior to downstream steps.")
    p.add_argument("--vcf-folder", dest="vcf_folder", default=None, help="Optional directory containing VCFs (overrides manifest paths).")
    return p.parse_args()


from contextlib import nullcontext  # at top of file

def main():
    args = parse_args()
    project_root = Path(args.folder).resolve()
    out_root = Path(args.out_dir).resolve()
    ensure_dir(out_root / "dna")

    # Read manifest
    df = read_table_guess(Path(args.manifest))

    # Optional simplified list
    if args.make_simplified:
        simp_out = Path(args.simplified) if args.simplified else (out_root / "case_ids.txt")
        outp = emit_simplified_case_list(Path(args.manifest), simp_out)
        logging.info(f"Wrote simplified Case-ID list: {outp}")

    # Normalize headers (case-insensitive lookup)
    cols = {c.lower().strip(): c for c in df.columns}

    case_col = None
    for key in ["case id", "case-id", "caseid"]:
        if key in cols:
            case_col = cols[key]
            break
    if case_col is None:
        raise SystemExit("dna_manifest.tsv must include a 'Case ID' column")

    file_id_col = None
    for key in ["file id", "file-id", "fileid"]:
        if key in cols:
            file_id_col = cols[key]
            break

    file_name_col = None
    for key in ["file name", "file-name", "filename"]:
        if key in cols:
            file_name_col = cols[key]
            break
    if file_name_col is None:
        raise SystemExit("dna_manifest.tsv must include a 'File Name' column")

    produced = 0
    seen = set()

    # Scratch context (only if --unzip-inputs enabled)
    scratch_mgr = _scratch(args) if getattr(args, "unzip_inputs", False) else nullcontext(Path())
    with scratch_mgr as scratch_dir:
        globals()['_SCRATCH_DIR'] = scratch_dir

        for _, row in df.iterrows():
            # Case-ID (strip anything after a comma)
            case_id = str(row[case_col]).split(",")[0].strip().strip('"')
            if case_id in seen:
                continue
            seen.add(case_id)
            safe_case = case_id.replace("/", "_")

            # Resolve file_id/file_name (file_id may be absent in some manifests)
            file_id = str(row[file_id_col]).strip() if file_id_col and pd.notna(row[file_id_col]) else ""
            file_name = str(row[file_name_col]).strip()
            if not file_name:
                logging.warning(f"[{case_id}] Missing 'File Name' in manifest; skipping")
                continue

            # Locate the source file under project_root
            try:
                src = resolve_dna_file(project_root, file_id, file_name)
            except FileNotFoundError as e:
                logging.warning(f"[{case_id}] {e}")
                continue

            # If requested, materialize compressed inputs to scratch
            eff_src = materialize_input(src, scratch_dir) if getattr(args, "unzip_inputs", False) else src

            # Parse VCF → DataFrame
            try:
                vdf = parse_vcf_to_df(eff_src, max_records=args.max_records)
            except Exception as e:
                logging.warning(f"[{case_id}] Failed to parse VCF '{eff_src}': {e}")
                continue

            # Write per-case CSV (filename = case only)
            vdf.insert(0, "Case-ID", case_id)
            out_path = out_root / "dna" / f"{safe_case}.csv"
            ensure_dir(out_path.parent)
            vdf.to_csv(out_path, index=False)
            produced += 1
            logging.info(f"[{case_id}] Wrote {out_path} ({len(vdf)} rows)")

    if produced == 0:
        logging.warning("No per-case DNA CSVs were produced. Check your dna_manifest.tsv and input files.")
    else:
        logging.info(f"Done. Wrote {produced} case CSVs under {out_root/'dna'}.")


if __name__ == "__main__":
    main()
