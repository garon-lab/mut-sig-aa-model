#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
DNA ANALYSIS

This script uses per-case DNA CSVs made from the DNA PREPROCESSOR step. It analyzes single
nucleotide variants (SNVs) and single nucleotide polymorphisms (SNPs), but can be adapted
as needed by the user. Expected input is MuTect-style flat TSVs with columns mapped to
['Case-ID', 'CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','NORMAL','TUMOR'].

It processes per-case DNA CSVs to:
1. Provide a summary of rows, SNVs and SNPs (or custom labels).
2. Create a signature composition (12-channel REF→ALT) summary across both groups.
3. Extract amino acid substitutions using the 16th INFO pipe-field (info.split('|')[15]) with
   first/last char as AA start/end, labeled by group.
4. Create 21×21 amino acid substitution matrices separately by group.

Notes:
- Default groups:
    SNP → FILTER contains "alt" and INFO contains "missense"
    SNV → FILTER contains "PASS" and INFO contains "missense"
- If you pass --labels X Y, those strings will replace "SNV"/"SNP" in outputs.

Usage:
python dna_analysis.py \
  --in_dir <input directory with dna/{case-id}.csv> \
  --simplified <case_ids.txt> \
  --out_dir <output directory> \
  [--summarize-variants --write-signatures --extract-mutations --write-matrices] \
  [--labels Foo Bar] [--max-records N] [--jobs N]
"""

from concurrent.futures import ProcessPoolExecutor, as_completed
import os
import argparse
import logging
from pathlib import Path
from typing import Iterable, Optional, List, Tuple
import pandas as pd

REQ_COLS = ['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','NORMAL','TUMOR']
AA_CODES = list("ACDEFGHIKLMNPQRSTVWY*")   # 21 AAs incl. stop
SUBS_12  = ["AC","AG","AT","CA","CG","CT","GA","GC","GT","TA","TC","TG"]

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

# ---------- helpers ----------
def _default_workers(jobs: Optional[int]) -> int:
    return jobs if (jobs and jobs > 0) else min(8, os.cpu_count() or 1)

def _locate_case_file(base_dir: Path, cid: str) -> Path:
    """
    Return the path to {cid}.csv whether base_dir is the DNA folder itself
    or its parent project folder.
    Tries: <base>/<cid>.csv, then <base>/dna/<cid>.csv
    """
    candidates = [
        base_dir / f"{cid}.csv",
        base_dir / "dna" / f"{cid}.csv",
    ]
    for p in candidates:
        if p.exists():
            return p
    # helpful error message listing both attempts
    raise FileNotFoundError(f"{candidates[0]} | {candidates[1]}")

def _dna_base(base_dir: Path) -> Path:
    """Return the directory that actually holds per-case DNA CSVs."""
    # Prefer when base is already the dna folder
    if (base_dir / "dna").is_dir():
        return base_dir / "dna"
    return base_dir

def ensure_dir(p: Path) -> None: p.mkdir(parents=True, exist_ok=True)
def load_case_ids(path: Path) -> List[str]: return [l.strip() for l in path.read_text().splitlines() if l.strip()]

def _read_case_csv(base_dir: Path, cid: str, max_records: Optional[int]) -> pd.DataFrame:
    """
    Read per-case DNA CSV (tolerates --in_dir being either the dna dir or its parent).
    Normalizes key columns:
      - unify 'Case-ID' header variants to 'Case-ID'
      - unify CHROM/POS/REF/ALT/FILTER/INFO casing
    """
    p = _locate_case_file(base_dir, cid)
    df = pd.read_csv(p, sep=None, engine="python", dtype=str).fillna("")

    # normalize 'Case-ID'
    lower_map = {c.lower().strip(): c for c in df.columns}
    case_key = next((lower_map[k] for k in ("case-id", "case id", "caseid") if k in lower_map), None)
    if case_key and case_key != "Case-ID":
        df = df.rename(columns={case_key: "Case-ID"})
        lower_map["case-id"] = "Case-ID"  # keep map in sync

    # optional cap
    if max_records is not None and len(df) > max_records:
        df = df.iloc[:max_records].copy()

    # require core fields (by lowercase)
    needed = ["chrom", "pos", "ref", "alt", "filter", "info"]
    missing = [k for k in needed if k not in lower_map]
    if missing:
        raise ValueError(f"{p} missing required columns: {', '.join(missing)}")

    # unify canonical casing for downstream access
    ren = {}
    for want, alts in {
        "CHROM": ["chrom"],
        "POS":   ["pos"],
        "REF":   ["ref"],
        "ALT":   ["alt"],
        "FILTER":["filter","filters"],
        "INFO":  ["info","ann","vep","csq","consequence","info.vep"],
    }.items():
        for k in alts:
            if k in lower_map:
                ren[lower_map[k]] = want
                break
    df = df.rename(columns=ren)

    return df


def _classify_label(df: pd.DataFrame, labels: Tuple[str,str]) -> pd.Series:
    f = df["FILTER"].astype(str); i = df["INFO"].astype(str)
    is_miss = i.str.contains("missense", case=False, na=False)
    lab = pd.Series("", index=df.index, dtype="string")
    lab[f.str.contains("alt", case=False, na=False) & is_miss]  = labels[1]  # default SNP
    lab[f.str.contains("PASS", case=False, na=False) & is_miss] = labels[0]  # default SNV
    return lab

def _first_alt_base(alt_col: pd.Series) -> pd.Series: return alt_col.astype(str).str.split(",").str[0].str.upper()
def _sbs_codes(ref: pd.Series, alt: pd.Series) -> pd.Series:
    mask = (ref.str.len()==1) & (alt.str.len()==1)
    return (ref[mask] + alt[mask]).str.replace(r"[^ACGT]", "", regex=True)

def _parse_aa_pair(info: str) -> Tuple[str,str]:
    try:
        field = str(info).split("|")[15]
        return (field[0], field[-1]) if field else ("","")
    except: return "",""

def _aa_matrix_from_pairs(pairs: Iterable[Tuple[str,str]]) -> pd.DataFrame:
    mat = pd.DataFrame(0, index=AA_CODES, columns=AA_CODES, dtype=int)
    for a0,a1 in pairs:
        if a0 in mat.index and a1 in mat.columns: mat.loc[a0,a1]+=1
    return mat

# ---------- 1) summary ----------
def _summary_one(base_dir: Path, cid: str, max_records: Optional[int], labels: Tuple[str,str]):
    try:
        df = _read_case_csv(base_dir, cid, max_records)
        f = df["FILTER"].astype(str); i = df["INFO"].astype(str)
        is_miss = i.str.contains("missense", case=False, na=False)
        is_lbl2 = f.str.contains("alt",  case=False, na=False) & is_miss   # SNP default
        is_lbl1 = f.str.contains("PASS", case=False, na=False) & is_miss   # SNV default
        return {"Case-ID": cid, labels[0]: int(is_lbl1.sum()), labels[1]: int(is_lbl2.sum()), "rows": int(len(df))}
    except Exception as e:
        logging.warning(f"[summary] {cid}: {e}")
        return None

def summarize_variants(in_dir: Path, case_ids: List[str], out_dir: Path,
                       max_records: Optional[int], labels: Tuple[str,str], jobs: Optional[int] = None):
    rows = []
    with ProcessPoolExecutor(max_workers=_default_workers(jobs)) as ex:
        futs = [ex.submit(_summary_one, in_dir, cid, max_records, labels) for cid in case_ids]
        for fut in as_completed(futs):
            r = fut.result()
            if r: rows.append(r)

    df_out = pd.DataFrame(rows)
    if "Case-ID" in df_out.columns:
        df_out = df_out.sort_values("Case-ID")
    outp = out_dir / "summary.csv"
    df_out.to_csv(outp, index=False)
    logging.info(f"[summary] wrote {outp}")


# ---------- 2) signatures ----------
def _sig_vec_from_df(df: pd.DataFrame) -> dict:
    """Build a 12-channel REF->ALT proportion vector from a variants DataFrame."""
    ref  = df["REF"].astype(str).str.upper()
    alt  = df["ALT"].astype(str).str.split(",").str[0].str.upper()
    mask = (ref.str.len()==1) & (alt.str.len()==1)
    sbs  = (ref[mask] + alt[mask]).str.replace(r"[^ACGT]", "", regex=True)
    counts = sbs.value_counts()
    vec = {k: int(counts.get(k, 0)) for k in SUBS_12}
    tot = max(1, sum(vec.values()))
    return {k: vec[k] / tot for k in SUBS_12}

def _signature_one(base_dir: Path, cid: str, max_records: Optional[int], labels: Tuple[str,str]):
    """
    Returns a tuple of three dicts (combined, label1, label2), each with Case-ID + 12-channel proportions.
    """
    try:
        df = _read_case_csv(base_dir, cid, max_records)

        # classify rows
        f = df["FILTER"].astype(str)
        i = df["INFO"].astype(str)
        is_miss = i.str.contains("missense", case=False, na=False)
        lab = pd.Series("", index=df.index, dtype="string")
        lab[f.str.contains("alt", case=False, na=False) & is_miss]  = labels[1]  # SNP
        lab[f.str.contains("PASS", case=False, na=False) & is_miss] = labels[0]  # SNV

        # subsets
        use_all = df[lab.isin(labels)]
        use_l1  = df[lab == labels[0]]
        use_l2  = df[lab == labels[1]]

        # build vectors (zeros if empty)
        row_all = {"Case-ID": cid, **({k:0.0 for k in SUBS_12} if use_all.empty else _sig_vec_from_df(use_all))}
        row_l1  = {"Case-ID": cid, **({k:0.0 for k in SUBS_12} if use_l1.empty  else _sig_vec_from_df(use_l1))}
        row_l2  = {"Case-ID": cid, **({k:0.0 for k in SUBS_12} if use_l2.empty  else _sig_vec_from_df(use_l2))}

        return row_all, row_l1, row_l2
    except Exception as e:
        logging.warning(f"[signatures] {cid}: {e}")
        zero = {"Case-ID": cid, **{k:0.0 for k in SUBS_12}}
        return zero, zero, zero

def write_signatures(in_dir: Path, case_ids: List[str], out_dir: Path, label_name: str,
                     max_records: Optional[int], labels: Tuple[str,str], jobs: Optional[int] = None):
    """
    Writes three files:
      <out_dir>/<LABEL1>-signature.csv
      <out_dir>/<LABEL2>-signature.csv
      <out_dir>/combined-signature.csv
    'label_name' is retained for compatibility but not used for file naming anymore.
    """
    rows_all, rows_l1, rows_l2 = [], [], []
    with ProcessPoolExecutor(max_workers=_default_workers(jobs)) as ex:
        futs = [ex.submit(_signature_one, in_dir, cid, max_records, labels) for cid in case_ids]
        for fut in as_completed(futs):
            r_all, r_l1, r_l2 = fut.result()
            rows_all.append(r_all); rows_l1.append(r_l1); rows_l2.append(r_l2)

    def _write(rows: List[dict], path: Path):
        df = pd.DataFrame(rows)
        if "Case-ID" in df.columns:
            df = df.sort_values("Case-ID")
        df.to_csv(path, index=False)

    # per-label files
    _write(rows_l1, out_dir / f"{labels[0]}-signature.csv")
    _write(rows_l2, out_dir / f"{labels[1]}-signature.csv")
    # combined file
    _write(rows_all, out_dir / "combined-signature.csv")

    logging.info(f"[signatures] wrote {out_dir / f'{labels[0]}-signature.csv'}, "
                 f"{out_dir / f'{labels[1]}-signature.csv'}, and {out_dir / 'combined-signature.csv'}")


# ---------- 3) amino acids ----------
def _aa_pair(info: str) -> Tuple[str,str]:
    try:
        fld = str(info).split("|")[15]
        return (fld[0] if fld else "", fld[-1] if fld else "")
    except Exception:
        return "",""

def _aa_one(base_dir: Path, cid: str, max_records: Optional[int], labels: Tuple[str,str]) -> Optional[pd.DataFrame]:
    try:
        df = _read_case_csv(base_dir, cid, max_records)
        f = df["FILTER"].astype(str); i = df["INFO"].astype(str)
        is_miss = i.str.contains("missense", case=False, na=False)
        lab = pd.Series("", index=df.index, dtype="string")
        lab[f.str.contains("alt", case=False, na=False) & is_miss]  = labels[1]
        lab[f.str.contains("PASS", case=False, na=False) & is_miss] = labels[0]
        use = df[lab.isin(labels)].copy()
        if use.empty:
            return None
        pairs = [ _aa_pair(s) for s in use["INFO"] ]
        use["AA_Start"] = [a for a,_ in pairs]
        use["AA_End"]   = [b for _,b in pairs]
        use["VarLabel"] = lab[use.index].values
        use["Case-ID"]  = cid
        return use[["Case-ID","VarLabel","CHROM","POS","REF","ALT","AA_Start","AA_End"]]
    except Exception as e:
        logging.warning(f"[aa] {cid}: {e}")
        return None

def extract_amino_acids(in_dir: Path, case_ids: List[str], out_dir: Path,
                        max_records: Optional[int], labels: Tuple[str,str], jobs: Optional[int] = None):
    parts: List[pd.DataFrame] = []
    with ProcessPoolExecutor(max_workers=_default_workers(jobs)) as ex:
        futs = [ex.submit(_aa_one, in_dir, cid, max_records, labels) for cid in case_ids]
        for fut in as_completed(futs):
            df = fut.result()
            if df is not None:
                parts.append(df)

    outp = out_dir / "aa-subs.csv"
    if parts:
        pd.concat(parts, ignore_index=True).to_csv(outp, index=False)
    else:
        pd.DataFrame(columns=["Case-ID","VarLabel","CHROM","POS","REF","ALT","AA_Start","AA_End"]).to_csv(outp, index=False)
    logging.info(f"[aa] wrote {outp}")



# ---------- 4) matrices ----------
def _matrix_one(base_dir: Path, cid: str, max_records: Optional[int],
                lbl: str, labels: Tuple[str,str], lbl_dir: Path) -> Optional[str]:
    try:
        df = _read_case_csv(base_dir, cid, max_records)
        lab = _classify_label(df, labels)
        use = df[lab == lbl]
        pairs = [_parse_aa_pair(s) for s in use["INFO"]]
        mat = _aa_matrix_from_pairs(pairs)
        outp = lbl_dir / f"{cid}.csv"
        mat.to_csv(outp)
        return str(outp)
    except Exception as e:
        logging.warning(f"[mat] {cid} ({lbl}): {e}")
        return None

def write_matrices(in_dir: Path, case_ids: List[str], out_dir: Path,
                   max_records: Optional[int], labels: Tuple[str,str], jobs: Optional[int] = None):
    workers = _default_workers(jobs)
    for lbl in labels:
        lbl_dir = out_dir / f"{lbl}-matrices"
        ensure_dir(lbl_dir)

        written = 0
        with ProcessPoolExecutor(max_workers=workers) as ex:
            futs = [ex.submit(_matrix_one, in_dir, cid, max_records, lbl, labels, lbl_dir) for cid in case_ids]
            for fut in as_completed(futs):
                r = fut.result()
                if r:
                    written += 1
        logging.info(f"[mat] wrote {written} matrices for label {lbl} -> {lbl_dir}")

# ---------- CLI ----------
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--in_dir",required=True); ap.add_argument("--simplified",required=True); ap.add_argument("--out_dir",required=True)
    ap.add_argument("--max-records",type=int,default=None); ap.add_argument("--jobs",type=int,default=None)
    ap.add_argument("--labels",nargs=2,metavar=("LABEL1","LABEL2"),help="Override default SNV/SNP labels")
    ap.add_argument("--summarize-variants",action="store_true"); ap.add_argument("--write-signatures",action="store_true")
    ap.add_argument("--extract-mutations",action="store_true"); ap.add_argument("--write-matrices",action="store_true")
    args=ap.parse_args()

    dna_dir_detected = _dna_base(Path(args.in_dir))
    logging.info(f"[analysis] Using DNA directory: {dna_dir_detected}")
    base = dna_dir_detected

    in_dir, out_dir = Path(args.in_dir), Path(args.out_dir); ensure_dir(out_dir)
    case_ids=load_case_ids(Path(args.simplified))
    labels=tuple(args.labels) if args.labels else ("SNV","SNP")

    if args.summarize_variants:
        summarize_variants(in_dir, case_ids, out_dir, args.max_records, labels, args.jobs)
    if args.write_signatures:
        write_signatures(in_dir, case_ids, out_dir, labels[0], args.max_records, labels, args.jobs)
    if args.extract_mutations:
        extract_amino_acids(in_dir, case_ids, out_dir, args.max_records, labels, args.jobs)
    if args.write_matrices:
        write_matrices(in_dir, case_ids, out_dir, args.max_records, labels, args.jobs)


if __name__=="__main__": main()
