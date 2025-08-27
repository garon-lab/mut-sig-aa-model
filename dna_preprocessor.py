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
  --in_dir <input directory> \
  --manifest <gdc-manifest tsv> \
  --out_dir <output directory>

(Optional)\
  --make-simplified \
  --preprocess-mutect \
  --vcf-folder <vcf directory> \
  --simplified <case_ids.txt> \
  --write-signatures \
  --extract-mutations \ 
  --write-matrices \
  --label <text>

Arguments:
Required
   --manifest           GDC-like TSV/CSV with at least Case ID, File Name (File ID optional)
   --in_dir             Input directory that contains raw data, format <folder>/dna/<File ID>/<File Name>
   --out_dir            Output directory that will contain <dna/<Case-ID>.csv that can be used in multiomic integration

General
   --max-records N      Cap parsed VCF rows per case (for smoke tests)
   --jobs               Controls parallel execution, if not provided, script uses min(8, CPU count)

Make/list Case-IDS:
   --make-simplified    Provides unique Case-IDs derived from --manifest
   --simplified         Path to write the Case-ID list (default: <out_dir>/case_ids.txt)
   
Preprocess: 
   --preprocess-mutect  Flag for extended analysis, strips '##' headers and writes prep/<Case-ID>.txt
   --vcf-folder         Where per-case VCFs live (default <folder/dna>)

Analytics (require --simplified file listing Case-IDs):\
   --simplified FILE              Path to case_ids.txt if it has been previously made, should have one Case-ID per line (no header)\
   --summarize-variants           Write SNP/SNV counts to <out.dir>/summary.csv\
   --write-signatures             Write <out_dir>/<label>-signature.csv\
   --label                        Label for file prefix (default: snv)\
   --extract-mutations            Extracts ST/END AA pairs to <out_dir>/<type>/<Case-ID>.csv\
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
import re

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

# Helpers

# Expected Mutect/VCF-like columns (no extra leading ID)
EXPECTED_MUTECT_COLS = ['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','NORMAL','TUMOR']

AA_LIST = ["A","R","N","D","C","Q","E","G","H","I",
           "L","K","M","F","P","S","T","W","Y","V"]
AA_SET = set(AA_LIST + ["*"])

AA3_TO_1 = {
    "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C",
    "GLN":"Q","GLU":"E","GLY":"G","HIS":"H","ILE":"I",
    "LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P",
    "SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V","TER":"*","STOP":"*"
}

def _is_single_base(s: str) -> bool:
    s = str(s).upper()
    return len(s) == 1 and s in {"A","C","G","T"}

def _normalize_one_letter(aa: str) -> str | None:
    if not isinstance(aa, str):
        return None
    s = aa.strip().upper()
    if s in AA_SET: 
        return s
    if len(s) == 3 and s in AA3_TO_1:
        return AA3_TO_1[s]
    return None

def _read_case_list(case_list_path: Path) -> list[str]:
    """Read a simple case list (handles header or raw IDs)."""
    ids = []
    for line in case_list_path.read_text().splitlines():
        s = line.strip().strip('"')
        if not s:
            continue
        # if it's TSV/CSV, take the first field
        if "\t" in s:
            s = s.split("\t", 1)[0]
        elif "," in s:
            s = s.split(",", 1)[0]
        # skip header-y lines
        if s.lower() in {"id", "case-id", "case id", "caseid"}:
            continue
        ids.append(s)
    return list(dict.fromkeys(ids))  # de-dupe, preserve order

def _make_unique_columns(cols):
    """Ensure duplicate column names become unique: ID, ID.1, ID.2, ..."""
    seen = {}
    out = []
    for c in cols:
        k = str(c)
        if k in seen:
            seen[k] += 1
            out.append(f"{k}.{seen[k]}")
        else:
            seen[k] = 0
            out.append(k)
    return out

def _read_mutect_or_vep_table(path: Path) -> pd.DataFrame:
    """
    Robust reader for prep files:
      - Headerless 11 cols  -> assign EXPECTED_MUTECT_COLS
      - Headerless 12 cols  -> ['EXTRA_ID'] + EXPECTED_MUTECT_COLS
      - Headered            -> use header; if still odd, take last 11 cols as the core set
    """
    # Try headerless first
    df = pd.read_csv(path, sep="\t", comment="#", header=None, engine="python")
    if df.empty:
        return df

    n = df.shape[1]
    if n == 11:
        df.columns = EXPECTED_MUTECT_COLS
    elif n == 12:
        # e.g., an extra first column with an external ID
        df.columns = ['EXTRA_ID'] + EXPECTED_MUTECT_COLS
    else:
        # Maybe there is a header row — try reading with header
        try:
            df2 = pd.read_csv(path, sep="\t", comment="#", header=0, engine="python")
            if set(EXPECTED_MUTECT_COLS).issubset(set(df2.columns)):
                df = df2
            else:
                # Last resort: take the last 11 cols as the core Mutect fields
                df = df.iloc[:, -11:]
                df.columns = EXPECTED_MUTECT_COLS
        except Exception:
            df = df.iloc[:, -11:]
            df.columns = EXPECTED_MUTECT_COLS

    # Avoid duplicate label issues (e.g., two "ID" columns)
    df.columns = _make_unique_columns(list(df.columns))
    return df

def _pick_col(df: pd.DataFrame, candidates):
    for c in candidates:
        if c in df.columns:
            return c
    return None

def _missense_mask(series: pd.Series) -> pd.Series:
    # Match both 'missense' and 'missense_variant' (VEP)
    return series.astype(str).str.contains(r"\bmissense\b|missense_variant", case=False, na=False)

def summarize_variants(case_list_path: Path, out_root: Path) -> Path:
    """
    Summarize SNP/SNV counts per case, preserving legacy logic:
      - SNP_alt_missense: FILTER contains 'alt'  AND INFO contains 'missense'
      - SNV_PASS_missense: FILTER contains 'PASS' AND INFO contains 'missense'
    Priority source: out_root/prep/{Case-ID}.txt (Mutect-like TSV)
    Fallback source: out_root/dna/{Case-ID}.csv (try same logic; else row count)
    Output: out_root/summary.csv with columns:
        Case-ID,SNP_alt_missense,SNV_PASS_missense,total_rows,source
    """
    prep_dir = out_root / "prep"
    dna_dir  = out_root / "dna"
    out_path = out_root / "summary.csv"

    if not case_list_path.exists():
        raise FileNotFoundError(f"Case list not found: {case_list_path}")

    case_ids = _read_case_list(case_list_path)
    if not case_ids:
        raise ValueError(f"No case IDs found in {case_list_path}")

    rows = []
    used_prep = used_dna = skipped = 0

    def _count_on_df(df: pd.DataFrame, source_label: str):
        total_rows = len(df)
        # Try to identify FILTER & INFO-like columns
        filter_col = _pick_col(df, ["FILTER","Filter","filter","FILTERS"])
        info_col   = _pick_col(df, ["INFO","Info","info","CSQ","Consequence","INFO.VEP","ANN","VEP"])
        if filter_col is not None and info_col is not None:
            filt = df[filter_col].astype(str)
            info = df[info_col].astype(str)
            snp_alt  = (filt.str.contains("alt",  case=False, na=False) & _missense_mask(info)).sum()
            snv_pass = (filt.str.contains("PASS", case=False, na=False) & _missense_mask(info)).sum()
            return int(snp_alt), int(snv_pass), int(total_rows), source_label
        # If missing, synthesize an INFO-like string from common annotation cols
        if info_col is None:
            cand = [c for c in df.columns if c.lower() in {"info","csq","consequence","ann","vep","vep_info"}]
            if cand:
                info = df[cand].astype(str).agg(" ".join, axis=1)
            else:
                info = pd.Series([""] * len(df))
        else:
            info = df[info_col].astype(str)
        if filter_col is None:
            filt = pd.Series([""] * len(df))
        else:
            filt = df[filter_col].astype(str)
        snp_alt  = (filt.str.contains("alt",  case=False, na=False) & _missense_mask(info)).sum()
        snv_pass = (filt.str.contains("PASS", case=False, na=False) & _missense_mask(info)).sum()
        return int(snp_alt), int(snv_pass), int(total_rows), source_label

    for cid in case_ids:
        safe = cid.replace("/", "_")
        prep_fp = prep_dir / f"{safe}.txt"
        dna_fp  = dna_dir  / f"{safe}.csv"

        # 1) Preferred: parse prep TSV (Mutect-like), tolerant of extra leading ID col
        if prep_fp.exists():
            try:
                vcf_df = _read_mutect_or_vep_table(prep_fp)
                snp_alt, snv_pass, total_rows, src = _count_on_df(vcf_df, "prep")
                rows.append({"Case-ID": cid,
                             "SNP_alt_missense": snp_alt,
                             "SNV_PASS_missense": snv_pass,
                             "total_rows": total_rows,
                             "source": src})
                used_prep += 1
                continue
            except Exception as e:
                logging.warning(f"[{cid}] Failed to parse prep file {prep_fp}: {e} (will try DNA fallback)")

        # 2) Fallback: per-case DNA CSV
        if dna_fp.exists():
            try:
                dfd = pd.read_csv(dna_fp)
                snp_alt, snv_pass, total_rows, src = _count_on_df(dfd, "dna")
                rows.append({"Case-ID": cid,
                             "SNP_alt_missense": snp_alt,
                             "SNV_PASS_missense": snv_pass,
                             "total_rows": total_rows,
                             "source": src})
                used_dna += 1
                continue
            except Exception as e:
                logging.warning(f"[{cid}] Failed to read {dna_fp}: {e}")

        logging.warning(f"Skipping {cid}: no usable prep or dna file under {out_root}")
        skipped += 1

    # Write output
    if not rows:
        logging.warning("No cases could be summarized. Writing empty summary.")
        pd.DataFrame(columns=["Case-ID","SNP_alt_missense","SNV_PASS_missense","total_rows","source"]).to_csv(out_path, index=False)
        return out_path

    out_df = pd.DataFrame(rows, columns=["Case-ID","SNP_alt_missense","SNV_PASS_missense","total_rows","source"])
    out_df.to_csv(out_path, index=False)
    logging.info(f"Summary written -> {out_path} (from prep: {used_prep}, from dna: {used_dna}, skipped: {skipped})")
    return out_path

def write_signatures(prep_dir: Path, simplified: Path, out_dir: Path, label: str):
    """
    Compute simple base-level mutation 'signatures' per case from preprocessed VCFs.
    Robust to:
      - headerless/headered prep files
      - optional extra leading 'ID' column
    Output: <out_dir>/<label>-signature.csv with header:
        Case-ID,SUM,CTGA,CAGT,GCCG,ATTA,AGTC,ACTG
    Logic:
      - Keep rows where FILTER contains 'alt' (case-insensitive)
      - SNVs only (REF/ALT are one of A/C/G/T and of length 1)
      - Count six folded categories:
           CTGA: CT + GA
           CAGT: CA + GT
           GCCG: CG + GC
           ATTA: AT + TA
           AGTC: AG + TC
           ACTG: AC + TG
    """
    prep_dir = Path(prep_dir); out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_file = out_dir / f"{label}-signature.csv"

    # Case list (tolerant to headers / TSV / CSV)
    sample_ids = _read_case_list(Path(simplified))
    header = ['Case-ID','SUM','CTGA','CAGT','GCCG','ATTA','AGTC','ACTG']
    rows = []

    for case_id in sample_ids:
        fp = prep_dir / f"{case_id}.txt"
        if not fp.exists():
            logging.warning(f"{case_id}: prep file not found at {fp}")
            continue
        try:
            df = _read_mutect_or_vep_table(fp)
            if df.empty:
                rows.append([case_id, 0, 0, 0, 0, 0, 0, 0]); continue

            # Require REF/ALT single nucleotide and FILTER contains 'alt'
            if 'FILTER' not in df.columns or 'REF' not in df.columns or 'ALT' not in df.columns:
                logging.warning(f"{case_id}: required columns missing in {fp.name}; skipping")
                continue
            snv_mask = df['REF'].map(_is_single_base) & df['ALT'].map(_is_single_base)
            alt_mask = df['FILTER'].astype(str).str.contains("alt", case=False, na=False)
            x = df.loc[snv_mask & alt_mask, ['REF','ALT']].astype(str).applymap(str.upper)

            # Count pairs
            def cnt(frm, to): return int(((x['REF']==frm) & (x['ALT']==to)).sum())
            counts = {
                'CTGA': cnt('C','T') + cnt('G','A'),
                'CAGT': cnt('C','A') + cnt('G','T'),
                'GCCG': cnt('C','G') + cnt('G','C'),
                'ATTA': cnt('A','T') + cnt('T','A'),
                'AGTC': cnt('A','G') + cnt('T','C'),
                'ACTG': cnt('A','C') + cnt('T','G'),
            }
            total = sum(counts.values())
            rows.append([case_id, total, counts['CTGA'], counts['CAGT'], counts['GCCG'],
                         counts['ATTA'], counts['AGTC'], counts['ACTG']])
        except Exception as e:
            logging.warning(f"Failed to process {case_id} signatures: {e}")

    pd.DataFrame(rows, columns=header).to_csv(out_file, index=False)
    logging.info(f"Signatures written -> {out_file}")


def extract_mutations(prep_dir: Path, out_dir: Path, simplified: Path, label: str):
    """
    Extract amino-acid pairs per mutation label ('snp' or 'snv') from preprocessed VCFs.

    FILTER rule (driven by label):
      - label == 'snp' -> rows where FILTER contains 'alt'
      - label == 'snv' -> rows where FILTER contains 'PASS'

    Output: <out_dir>/<label>/<Case-ID>-<label>.csv with columns: ST, END, #CHROM, TUMOR
    """
    mtype = str(label).lower().strip()
    if mtype not in {"snp", "snv"}:
        logging.warning(f"extract_mutations: label '{label}' not in {{snp,snv}}; defaulting to 'snv'")
        mtype = "snv"

    prep_dir = Path(prep_dir); out_dir = Path(out_dir); simplified = Path(simplified)
    out_sub = out_dir / mtype
    out_sub.mkdir(parents=True, exist_ok=True)

    ids = _read_case_list(simplified)

    aa_pair_regex = re.compile(r"\b([ACDEFGHIKLMNPQRSTVWY\*])\s*/\s*([ACDEFGHIKLMNPQRSTVWY\*])\b")
    hgvsp_1L = re.compile(r"p\.([ACDEFGHIKLMNPQRSTVWY\*])\w*\d+([ACDEFGHIKLMNPQRSTVWY\*])", re.IGNORECASE)
    hgvsp_3L = re.compile(r"p\.([A-Z][a-z]{2})\d+([A-Z][a-z]{2})")

    def parse_info_to_pair(info: str):
        s = str(info)
        m = aa_pair_regex.search(s)
        if m:
            a, b = _normalize_one_letter(m.group(1)), _normalize_one_letter(m.group(2))
            if a in AA_SET and b in AA_SET: return a, b
        m = hgvsp_1L.search(s)
        if m:
            a, b = _normalize_one_letter(m.group(1)), _normalize_one_letter(m.group(2))
            if a in AA_SET and b in AA_SET: return a, b
        m = hgvsp_3L.search(s)
        if m:
            a, b = _normalize_one_letter(AA3_TO_1.get(m.group(1).upper())), _normalize_one_letter(AA3_TO_1.get(m.group(2).upper()))
            if a in AA_SET and b in AA_SET: return a, b
        return None

    for case_id in ids:
        try:
            infile = prep_dir / f"{case_id}.txt"
            if not infile.exists():
                logging.warning(f"extract_mutations: {case_id} missing prep file {infile}")
                continue

            vcf_df = _read_mutect_or_vep_table(infile)
            if vcf_df.empty:
                logging.warning(f"extract_mutations: {case_id} empty file {infile}")
                continue

            if 'FILTER' not in vcf_df.columns:
                logging.warning(f"{case_id}: no FILTER column in {infile.name}; skipping")
                continue

            filt_mask = (vcf_df['FILTER'].astype(str).str.contains("alt", case=False, na=False)
                         if mtype == "snp"
                         else vcf_df['FILTER'].astype(str).str.contains("PASS", case=False, na=False))

            sub = vcf_df.loc[filt_mask].copy()
            if sub.empty:
                logging.info(f"{case_id}: no rows after {mtype} filter")
                continue

            info_col = next((c for c in ["INFO","Info","info","CSQ","ANN","VEP","INFO.VEP"] if c in sub.columns), None)
            if info_col is None:
                logging.warning(f"{case_id}: no INFO/CSQ/ANN column; skipping AA extraction")
                continue

            pairs = sub[info_col].apply(parse_info_to_pair)
            keep = pairs.notna()
            if not keep.any():
                logging.warning(f"{case_id}: no parseable AA changes in INFO; skipping")
                continue

            st_end = pairs[keep].tolist()
            st = [a for a,_ in st_end]
            ed = [b for _,b in st_end]
            chrom = sub.loc[keep, '#CHROM'] if '#CHROM' in sub.columns else pd.Series([""]*len(st))
            tumor = sub.loc[keep, 'TUMOR'] if 'TUMOR' in sub.columns else pd.Series([""]*len(st))

            aa_df = pd.DataFrame({'ST': st, 'END': ed, '#CHROM': chrom.values, 'TUMOR': tumor.values})
            outfile = out_sub / f"{case_id}-{mtype}.csv"
            aa_df.to_csv(outfile, index=False)
            logging.info(f"[{case_id}] extracted {len(aa_df)} AA pairs -> {outfile}")
        except Exception as e:
            logging.warning(f"Failed to extract {mtype} for {case_id}: {e}")


def generate_aa_matrix(df: pd.DataFrame) -> pd.DataFrame:
    """
    Build a 21x21 matrix from ST/END columns (tolerant of 3-letter codes and lowercase).
    Missing rows/cols filled with 0. STOP is '*'.
    """
    if not {'ST','END'}.issubset(df.columns):
        raise ValueError("generate_aa_matrix: input DataFrame must have ST and END columns")

    # Normalize to one-letter uppercase (map 3-letter if present)
    st = df['ST'].apply(_normalize_one_letter)
    en = df['END'].apply(_normalize_one_letter)
    ok = st.notna() & en.notna()
    if not ok.any():
        # return a zero matrix with proper labels
        idx = AA_LIST + ["*"]; cols = AA_LIST + ["*"]
        return pd.DataFrame(0, index=idx, columns=cols, dtype=int)

    g = pd.DataFrame({'ST': st[ok], 'END': en[ok]})
    mat = g.groupby(['ST','END']).size().unstack(fill_value=0)
    mat = mat.reindex(index=AA_LIST+["*"], columns=AA_LIST+["*"], fill_value=0)
    return mat


def write_matrices(out_dir: Path, simplified: Path):
    """
    Write per-case amino-acid 21x21 matrices under:
        <out_dir>/snp/matrices/<Case-ID>.csv
        <out_dir>/snv/matrices/<Case-ID>.csv
    Expects extract_mutations() outputs at:
        <out_dir>/snp/<Case-ID>-snp.csv  and  <out_dir>/snv/<Case-ID>-snv.csv
    """
    out_dir = Path(out_dir)
    ids = _read_case_list(Path(simplified))
    for mtype in ['snp','snv']:
        (out_dir / mtype / "matrices").mkdir(parents=True, exist_ok=True)

    for case_id in ids:
        for mtype in ['snp','snv']:
            try:
                infile = out_dir / mtype / f"{case_id}-{mtype}.csv"
                if not infile.exists():
                    logging.info(f"write_matrices: missing {infile}, skipping")
                    continue
                df_aa = pd.read_csv(infile)
                mat = generate_aa_matrix(df_aa)
                out_path = out_dir / mtype / "matrices" / f"{case_id}.csv"
                mat.to_csv(out_path)
                logging.info(f"[{case_id}] {mtype} matrix -> {out_path}")
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
            pick_col = 0 if df.shape[1] > 0 else None
            if df.shape[1] > 5:
                pick_col = 5
            ids = df.iloc[:, pick_col].astype(str)
    except Exception:
        # Plain text list
        ids = [line.strip() for line in Path(manifest_file).read_text().splitlines() if line.strip()]

    for case_id in ids:
        case_id = str(case_id).split(",")[0].strip().strip('"')
        vcf_path = _find_vcf_for_case(folder, case_id)
        target = Path(out_dir) / "prep" / f"{case_id}.txt"
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
    p = argparse.ArgumentParser(
        description="Create per-case DNA CSVs (dna/{Case-ID}.csv) from a dna_manifest.tsv using a VCF parser."
    )
    p.add_argument("--folder", "--project-root", "--data-root", "--input_dir", "--in_dir", required=True,
                   help="Project root containing dna/ or vep/ subfolders")
    p.add_argument("--manifest", required=True,
                   help="Path to dna_manifest.tsv (GDC-like)")
    p.add_argument("--out_dir", required=True,
                   help="Output directory root (CSV written under out_dir/dna/)")
    p.add_argument("--unzip-inputs", "--unzip_inputs", dest="unzip_inputs",
                   action="store_true",
                   help="Automatically materialize compressed inputs (.gz, .zip) into a scratch directory before processing.")
    p.add_argument("--scratch-dir", "--scratch_dir", dest="scratch_dir",
                   default=None,
                   help="Optional directory to use for temporary extracted files. Defaults to a new Temporary Directory each run.")
    p.add_argument("--max-records", "--max_records", dest="max_records",
                   type=int, default=None,
                   help="Optional cap on parsed VCF records per case (for testing)")
    p.add_argument("--make-simplified", "--make_simplified", dest="make_simplified",
                   action="store_true",
                   help="Emit a Case-ID list derived from --manifest")
    p.add_argument("--simplified",
                   help="Path to write the Case-ID list (default: <out_dir>/case_ids.txt)")
    p.add_argument("--preprocess-mutect", "--preprocess_mutect", dest="preprocess_mutect",
                   action="store_true",
                   help="Preprocess Mutect VCFs prior to downstream steps.")
    p.add_argument("--vcf-folder", "--vcf_folder", dest="vcf_folder",
                   default=None,
                   help="Optional directory containing VCFs (overrides manifest paths).")
    p.add_argument("--jobs", type=int, default=None,
                   help="Parallel workers (currently not used; accepted for compatibility)")

    # Analytics controls
    p.add_argument("--summarize-variants", "--summarize_variants", dest="summarize_variants",
                   action="store_true",
                   help="Summarize SNP/SNV counts from out_dir/prep/*.txt into out_dir/summary.csv")
    p.add_argument("--write-signatures", "--write_signatures", dest="write_signatures",
                   action="store_true",
                   help="Write <out_dir>/<label>-signature.csv")
    p.add_argument("--extract-mutations", "--extract_mutations", dest="extract_mutations",
                   action="store_true",
                   help="Extract ST/END AA pairs to <out_dir>/<label>/<Case-ID>-<label>.csv "
                        "(uses label to choose SNP/SNV logic)")
    p.add_argument("--write-matrices", "--write_matrices", dest="write_matrices",
                   action="store_true",
                   help="Write 21x21 AA matrices to <out_dir>/{snp,snv}/matrices/<Case-ID>.csv")

    # Single knob that also drives mutation type
    p.add_argument("--label", dest="label",
                   default="snv",
                   help="Label used for filenames AND mutation-type logic in extraction "
                        "(accepted: 'snv' or 'snp'; default: snv)")
    return p.parse_args()



def main():
    args = parse_args()
    project_root = Path(args.folder).resolve()
    out_root = Path(args.out_dir).resolve()
    ensure_dir(out_root / "dna")

    # Read manifest
    df = read_table_guess(Path(args.manifest))

    # Path for the simplified list (compute once)
    simp_out = Path(args.simplified) if args.simplified else (out_root / "case_ids.txt")

    # Optional: emit simplified Case-ID list
    if args.make_simplified:
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

            # Resolve file_id/file_name (file_id may be absent)
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
          
            # Keep only PASS variants
            if "FILTER" in vdf.columns:
                vdf_pass = vdf[vdf["FILTER"] == "PASS"].copy()
            else:
                # Be robust if header differs in case
                cols_norm = {c.lower(): c for c in vdf.columns}
                if "filter" in cols_norm:
                    vdf_pass = vdf[vdf[cols_norm["filter"]] == "PASS"].copy()
                else:
                    logging.warning(f"[DNA] {case_id}: no FILTER column found; skipping write.")
                    vdf_pass = vdf.iloc[0:0].copy()  # empty
            
            # Only write if there’s at least one PASS variant
            if vdf_pass.empty:
                logging.warning(f"[DNA] {case_id}: no PASS variants — skipping dna/{case_id}.csv")
            else:
                vdf_pass.to_csv(out_path, index=False)
            

            vdf.to_csv(out_path, index=False)
            produced += 1
            logging.info(f"[{case_id}] Wrote {out_path} ({len(vdf)} rows)")

    if produced == 0:
        logging.warning("No per-case DNA CSVs were produced. Check your dna_manifest.tsv and input files.")
    else:
        logging.info(f"Done. Wrote {produced} case CSVs under {out_root/'dna'}.")

        # ---- Optional downstream steps ----
    if args.preprocess_mutect and 'preprocess_mutect' in globals():
        try:
            vcf_root = Path(args.vcf_folder) if args.vcf_folder else project_root
            preprocess_mutect(vcf_root, Path(args.manifest), out_root)
        except Exception as e:
            logging.error(f"preprocess_mutect failed: {e}")

    # ensure we have a case list for downstream steps
    if (args.summarize_variants or args.write_signatures or
        args.extract_mutations or args.write_matrices):
        if not simp_out.exists():
            try:
                emit_simplified_case_list(Path(args.manifest), simp_out)
                logging.info(f"Wrote simplified Case-ID list: {simp_out}")
            except Exception as e:
                logging.error(f"failed to produce case list for downstream steps: {e}")

    if args.summarize_variants and 'summarize_variants' in globals():
        try:
            summarize_variants(simp_out, out_root)
        except Exception as e:
            logging.error(f"summarize_variants failed: {e}")

    if args.write_signatures and 'write_signatures' in globals():
        try:
            write_signatures(out_root / "prep", simp_out, out_root, label=args.label)
        except Exception as e:
            logging.error(f"write_signatures failed: {e}")

    if args.extract_mutations and 'extract_mutations' in globals():
        try:
            extract_mutations(out_root / "prep", out_root, simp_out, label=args.label)
        except Exception as e:
            logging.error(f"extract_mutations failed: {e}")

    if args.write_matrices and 'write_matrices' in globals():
        try:
            write_matrices(out_root, simp_out)
        except Exception as e:
            logging.error(f"write_matrices failed: {e}")



if __name__ == "__main__":
    main()
