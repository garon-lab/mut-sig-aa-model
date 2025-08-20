#!/usr/bin/env python3
"""
Comparison Pipeline

Summarizes observed amino acid variant counts per sample and generates similarity
and visualization outputs.

Features
--------
1) Aggregate observed amino acid variants into one vector per sample (summary.csv)
   - Auto-detects 'tall' (ID, ST, END, Count) OR 21×21 matrix files in --observed-dir
2) Compare observed vectors vs. comparison/reference vectors (cosine similarity)
3) Visualize similarity as a heatmap
4) Plot clustermaps (counts & row-wise proportions)
5) Single-file utilities: plot counts/proportions and optional compare

Dependencies: pandas, numpy, matplotlib, seaborn

Typical usage:
    python comparison_and_modeling.py \
        --obs_dir <observed CSVs> \
        --comp_dir <directory containing comparison_vectors.csv> \
        --manifest <optional manifest (first col are IDs)> \
        --out_dir <output> \
        --step all

Single-file usage (you already have a summary CSV):
    python comparison_and_modeling.py \
        --vector_file path/to/summary.csv \
        --comp_dir <directory with comparison_vectors.csv> \
        --out_dir <output> \
        --step single-file
"""

import argparse
import logging
from pathlib import Path
from typing import Optional, Set, Dict
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sb

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

# ---------- Helpers & Constants ----------

AA3_TO_1: Dict[str, str] = {
    "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C",
    "GLN":"Q","GLU":"E","GLY":"G","HIS":"H","ILE":"I",
    "LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P",
    "SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V",
    "TER":"*", "STOP":"*"
}
AA1_SET = set(list("ARNDCQEGHILKMFPSTWYV") + ["X"])  # 'X' used for STOP here

def _normalize_id_col(df: pd.DataFrame) -> pd.DataFrame:
    rename_map = {
        'Case-ID': 'ID', 'case-id': 'ID', 'CaseId': 'ID',
        'CaseID': 'ID', 'caseID': 'ID', 'Case_Id': 'ID'
    }
    df = df.rename(columns=rename_map)
    if 'ID' not in df.columns:
        raise KeyError("Observed file is missing an 'ID' (or Case-ID) column.")
    return df

def normalize_aa_cols(df: pd.DataFrame) -> pd.DataFrame:
    """Normalize any reasonable AA start/end columns to ST/END."""
    candidates = [
        ('st','end'),
        ('fromaa','toaa'),
        ('from','to'),
        ('refaa','altaa'),
        ('ref_aa','alt_aa'),
    ]
    lower = {c.lower(): c for c in df.columns}
    for left, right in candidates:
        if left in lower and right in lower:
            return df.rename(columns={lower[left]:'ST', lower[right]:'END'})
    raise KeyError("Could not find AA start/end columns; expected one of "
                   "ST/END, FromAA/ToAA, from/to, RefAA/AltAA, Ref_AA/Alt_AA")

def _norm_symbol_to_one_letter_or_X(x: object) -> Optional[str]:
    """Map row/col headers to 1-letter AA; STOP/*/TER/Ter/X => 'X'."""
    if x is None:
        return None
    s = str(x).strip()
    if not s:
        return None
    su = s.upper()
    if su in ("*", "STOP", "TER", "TERM", "X"):
        return "X"
    if len(su) == 1 and su in AA1_SET.union(set("ARNDCQEGHILKMFPSTWYV")):
        # If plain 1-letter but not X (stop), return as-is; X already handled
        return su if su != "*" else "X"
    if len(su) == 3 and su in AA3_TO_1:
        one = AA3_TO_1[su]
        return "X" if one == "*" else one
    return None

# AA-pair column order for vectors
HEADING = ['AE','AG','AP','AS','AT','AV','CF','CG','CR','CS','CW','CX','CY',
           'DA','DE','DG','DH','DN','DV','DY','EA','ED','EG','EK','EQ','EV','EX','FC',
           'FI','FL','FS','FV','FY','GA','GC','GD','GE','GR','GS','GV','GW','GX','HD',
           'HL','HN','HP','HQ','HR','HY','IF','IK','IL','IM','IN','IR','IS','IT','IV',
           'IW','KE','KI','KM','KN','KQ','KR','KT','KX','LF','LH','LI','LM','LP','LQ',
           'LR','LS','LV','LW','LX','MI','MK','ML','MR','MT','MV','ND','NH','NI','NK',
           'NS','NT','NY','PA','PH','PL','PQ','PR','PS','PT','QE','QH','QK','QL','QP',
           'QR','QX','RC','RG','RH','RI','RK','RL','RM','RP','RQ','RS','RT','RW','RX',
           'SA','SC','SF','SG','SI','SL','SN','SP','SR','ST','SW','SX','SY','TA','TI',
           'TK','TM','TN','TP','TR','VA','VD','VE','VF','VG','VI','VL','VM','WC','WG',
           'WL','WR','WS','WX','XC','XE','XG','XK','XL','XQ','XR','XS','XW','XY','YC',
           'YD','YF','YH','YN','YS','YX']

# ---------- IO & Vectorization ----------

def _read_manifest_ids(manifest_path: Optional[str]) -> Optional[Set[str]]:
    if not manifest_path:
        return None
    m = pd.read_table(manifest_path, header=None).iloc[:, 0].astype(str)
    return set(m)

def _collect_observed_tall(observed_dir: str, manifest_ids: Optional[Set[str]]) -> pd.DataFrame:
    """
    Ingest all 'tall' observed CSVs: must contain ID, ST, END and (optionally) Count.
    If Count is missing, assume Count=1 for each row.
    """
    parts = []
    for fp in Path(observed_dir).glob("*.csv"):
        try:
            df = pd.read_csv(fp)
            df = _normalize_id_col(df)
            df = normalize_aa_cols(df)
            if 'Count' not in df.columns:
                df['Count'] = 1
            if manifest_ids is not None:
                df = df[df['ID'].astype(str).isin(manifest_ids)]
            parts.append(df[['ID','ST','END','Count']])
        except Exception as e:
            logging.debug(f"[tall adapter] skipping {fp.name}: {e}")
    if not parts:
        raise FileNotFoundError("No tall-format files found.")
    return pd.concat(parts, ignore_index=True)

def _vectorize_one_matrix_file(fp: Path) -> Optional[pd.DataFrame]:
    """
    Read a single 21×21-like matrix and return a 1-row DataFrame with columns in HEADING.
    File name stem is used as sample ID.
    """
    try:
        raw = pd.read_csv(fp, index_col=0)
    except Exception as e:
        logging.debug(f"[matrix adapter] cannot read {fp.name}: {e}")
        return None

    # Normalize axis labels and choose orientation that yields more valid AA labels
    rows_norm = [ _norm_symbol_to_one_letter_or_X(x) for x in raw.index ]
    cols_norm = [ _norm_symbol_to_one_letter_or_X(x) for x in raw.columns ]
    row_valid = sum(1 for r in rows_norm if r in AA1_SET)
    col_valid = sum(1 for c in cols_norm if c in AA1_SET)

    mat = raw.copy()
    if row_valid < col_valid:
        mat = mat.T
        rows_norm, cols_norm = cols_norm, rows_norm
        row_valid, col_valid = col_valid, row_valid

    # Apply normalized labels
    mat.index = rows_norm
    mat.columns = cols_norm
    # Keep only valid AAs
    mat = mat.loc[[r for r in mat.index if r in AA1_SET],
                  [c for c in mat.columns if c in AA1_SET]]
    if mat.empty:
        return None

    # Force numeric
    mat = mat.apply(pd.to_numeric, errors='coerce').fillna(0)

    # Build row for this sample
    sample_id = fp.stem
    row = {'ID': sample_id}
    for code in HEADING:
        st, en = code[0], code[1]
        val = 0
        if st in mat.index and en in mat.columns:
            v = mat.at[st, en]
            try:
                val = float(v)
            except Exception:
                val = 0
        row[code] = val
    row['SUM'] = sum(row[c] for c in HEADING)
    return pd.DataFrame([row], columns=['ID','SUM'] + HEADING)

def _collect_observed_matrices(observed_dir: str, manifest_ids: Optional[Set[str]]) -> pd.DataFrame:
    rows = []
    for fp in Path(observed_dir).glob("*.csv"):
        if manifest_ids and fp.stem not in manifest_ids:
            continue
        r = _vectorize_one_matrix_file(fp)
        if r is not None:
            rows.append(r)
    if not rows:
        raise FileNotFoundError("No matrix-format files could be vectorized.")
    return pd.concat(rows, ignore_index=True)

def summarize_observed(observed_dir: str, out_dir: str, manifest_path: Optional[str]) -> pd.DataFrame:
    """
    Produce out_dir/summary.csv
    Columns: ID, SUM, <AA-pairs in HEADING...>
    Auto-detect tall vs matrix inputs.
    """
    manifest_ids = _read_manifest_ids(manifest_path)

    # Try tall first, fall back to matrix adapter
    try:
        tall = _collect_observed_tall(observed_dir, manifest_ids)
        tall['code'] = tall['ST'].astype(str) + tall['END'].astype(str)
        piv = tall.pivot_table(index='ID', columns='code', values='Count',
                               aggfunc='sum', fill_value=0)
        for c in HEADING:
            if c not in piv.columns:
                piv[c] = 0
        piv = piv[HEADING]
        piv.insert(0, 'SUM', piv.sum(axis=1))
        out_df = piv.reset_index()
        logging.info("Summarized observed (tall format detected).")
    except Exception:
        # Matrix adapter
        out_df = _collect_observed_matrices(observed_dir, manifest_ids)
        # Ensure ordering
        out_df = out_df[['ID','SUM'] + HEADING]
        logging.info("Summarized observed (matrix format detected).")

    out_path = Path(out_dir) / "summary.csv"
    out_df.to_csv(out_path, index=False)
    logging.info(f"Saved summary -> {out_path}")
    return out_df

def load_comparison(comparison_dir: str, comparison_csv: Optional[str] = None) -> pd.DataFrame:
    """
    Load comparison/reference vectors. Default file name: comp_dir/comparison_vectors.csv
    Must have an ID (first column if unnamed) and AA-pair columns.
    Returns dataframe indexed by ID with only AA-pair columns (ordered by HEADING when present).
    """
    path = Path(comparison_csv) if comparison_csv else Path(comparison_dir) / "comparison_vectors.csv"
    df = pd.read_csv(path)
    if 'ID' not in df.columns:
        df = df.rename(columns={df.columns[0]: 'ID'})
    df = df.set_index('ID')

    cols_present = [c for c in HEADING if c in df.columns]
    if not cols_present:
        raise ValueError(f"Comparison vectors file '{path}' has no AA-pair columns matching the 21×21 codes.")
    if len(cols_present) < len(HEADING):
        logging.warning(f"Comparison vectors missing {len(HEADING) - len(cols_present)} AA codes (will align to overlap).")

    return df[cols_present]

# ---------- Comparisons & Plots ----------

def compare_vectors(observed_summary: pd.DataFrame, comparison_df: pd.DataFrame, out_dir: str) -> pd.DataFrame:
    """
    Cosine similarity between observed (rows=samples) and comparison (rows=reference labels).
    observed_summary must have columns: ID, SUM, <AA-pairs...>
    """
    obs = observed_summary.set_index('ID').drop(columns=['SUM'], errors='ignore')
    common = [c for c in obs.columns if c in comparison_df.columns]
    if not common:
        raise ValueError("No overlapping AA-pair columns between observed and comparison.")

    X = obs[common].to_numpy(dtype=float)
    Y = comparison_df[common].to_numpy(dtype=float)

    Xn = X / (np.linalg.norm(X, axis=1, keepdims=True) + 1e-12)
    Yn = Y / (np.linalg.norm(Y, axis=1, keepdims=True) + 1e-12)
    sims = Xn @ Yn.T

    sim_df = pd.DataFrame(sims, index=obs.index, columns=comparison_df.index)
    out_path = Path(out_dir) / "similarity_matrix.csv"
    sim_df.to_csv(out_path)
    logging.info(f"Saved similarity matrix -> {out_path}")
    return sim_df

def plot_similarity_heatmap(sim_df: pd.DataFrame, out_dir: str, name: str = "heatmap.png"):
    plt.figure(figsize=(10, 10))
    sb.heatmap(sim_df.astype(float), cmap='viridis')
    out_path = Path(out_dir) / name
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    plt.close()
    logging.info(f"Saved heatmap -> {out_path}")

def plot_clustermap(summary_df: pd.DataFrame, out_dir: str, proportions: bool = False, name: Optional[str] = None):
    """
    Clustermap of sample-by-AApair matrix (counts or row-wise proportions).
    summary_df columns: ID, SUM, AA-pairs...
    """
    mat = summary_df.set_index('ID').drop(columns=['SUM'], errors='ignore').copy()
    if proportions:
        denom = summary_df.set_index('ID')['SUM'].replace(0, np.nan)
        mat = mat.div(denom, axis=0).fillna(0) * 100.0  # percent
    cm = sb.clustermap(mat.T, metric='cosine', method='average', cmap='viridis',
                       figsize=(18, 22))
    out_path = Path(out_dir) / (name or ("clustermap_prop.png" if proportions else "clustermap_counts.png"))
    cm.savefig(out_path, dpi=200)
    plt.close(cm.fig)
    logging.info(f"Saved clustermap -> {out_path}")

# ---------- Single-file utilities ----------

def _read_summary_csv(vector_file: str) -> pd.DataFrame:
    """
    Expect a CSV with columns: ID, SUM, <AA-pair columns...>
    If SUM is missing, we compute it.
    """
    df = pd.read_csv(vector_file)
    if 'ID' not in df.columns:
        df = df.rename(columns={df.columns[0]: 'ID'})
    cols_present = [c for c in HEADING if c in df.columns]
    if 'SUM' not in df.columns:
        df['SUM'] = df[cols_present].sum(axis=1)
    # reorder
    return df[['ID','SUM'] + cols_present]

def single_file_plots(vector_file: str, out_dir: str):
    """
    From a precomputed vector CSV -> heatmaps and clustermaps (counts + proportions).
    """
    df = _read_summary_csv(vector_file)

    # counts heatmap
    mat_counts = df.set_index('ID').drop(columns=['SUM'], errors='ignore').T
    plt.figure(figsize=(18, 12))
    sb.heatmap(mat_counts, cmap='viridis')
    out_path = Path(out_dir) / "aa-counts.png"
    plt.tight_layout(); plt.savefig(out_path, dpi=200); plt.close()
    logging.info(f"Saved -> {out_path}")

    # proportions heatmap
    denom = df.set_index('ID')['SUM'].replace(0, np.nan)
    mat_prop = df.set_index('ID').drop(columns=['SUM'], errors='ignore').div(denom, axis=0).fillna(0).T * 100.0
    plt.figure(figsize=(18, 12))
    sb.heatmap(mat_prop, cmap='viridis', vmin=0, vmax=5)
    out_path = Path(out_dir) / "aa-proportions.png"
    plt.tight_layout(); plt.savefig(out_path, dpi=200); plt.close()
    logging.info(f"Saved -> {out_path}")

    # clustermaps
    plot_clustermap(df, out_dir, proportions=False, name="aa-clustermap-counts.png")
    plot_clustermap(df, out_dir, proportions=True,  name="aa-clustermap-proportions.png")

def single_file_compare(vector_file: str, comparison_dir: str, out_dir: str, comparison_csv: Optional[str] = None):
    """
    Compare vectors from a precomputed summary file to comparison vectors.
    Saves single_similarity.csv and single_similarity_heatmap.png
    """
    df = _read_summary_csv(vector_file)
    comp = load_comparison(comparison_dir, comparison_csv=comparison_csv)
    sim = compare_vectors(df, comp, out_dir)
    plot_similarity_heatmap(sim, out_dir, name="single_similarity_heatmap.png")

# ---------- CLI ----------

def parse_args():
    p = argparse.ArgumentParser(description="Comparison Pipeline")
    p.add_argument('--observed-dir', '--observed_dir', '--obs-dir', '--obs_dir',
                   dest='observed_dir',
                   help='Directory of observed AA CSVs (tall files or 21x21 matrices)')
    p.add_argument('--comparison-dir', '--comparison_dir', '--comp-dir', '--comp_dir',
                   dest='comparison_dir', required=False,
                   help='Directory containing comparison_vectors.csv (for compare/heatmap/single-file compare)')
    p.add_argument('--comparison-csv', dest='comparison_csv', default=None,
                   help="Optional explicit path to comparison/reference CSV (overrides default comparison_vectors.csv)")
    p.add_argument('--manifest', help='Optional manifest (first column are sample IDs) to filter observed')
    p.add_argument('--out-dir', '--out_dir', dest='out_dir', required=True,
                   help='Directory for outputs')
    p.add_argument('--vector-file', '--vector_file', dest='vector_file',
                   help='Path to an existing summary CSV (used by --step single-file or to skip summarize)')
    p.add_argument('--step', choices=['summarize','compare','heatmap','single-file','all'],
                   default='all', help='Which step(s) to run')
    return p.parse_args()

# ---------- Main ----------

def main():
    args = parse_args()
    out_dir = Path(args.out_dir); out_dir.mkdir(parents=True, exist_ok=True)

    summary_df = None
    sim_df = None

    # SUMMARIZE
    if args.step in ('summarize', 'all'):
        if not args.observed_dir:
            raise SystemExit("--observed-dir is required for step 'summarize' or 'all'")
        summary_df = summarize_observed(args.observed_dir, args.out_dir, args.manifest)

    # If skipping summarize but we need the summary later, load it
    if summary_df is None and args.step in ('compare','heatmap','all'):
        vf = args.vector_file or (Path(args.out_dir) / 'summary.csv')
        if not Path(vf).exists():
            raise SystemExit(f"Summary CSV not found: {vf}. Either run --step summarize first or pass --vector_file.")
        summary_df = _read_summary_csv(str(vf))

    # COMPARE
    if args.step in ('compare','all'):
        if not args.comparison_dir and not args.comparison_csv:
            raise SystemExit("--comparison-dir or --comparison-csv is required for 'compare'/'all'")
        comp = load_comparison(args.comparison_dir or "", comparison_csv=args.comparison_csv)
        sim_df = compare_vectors(summary_df, comp, args.out_dir)

    # HEATMAP
    if args.step in ('heatmap','all'):
        if sim_df is None:
            # load previously saved similarity if not computed in this run
            sim_path = Path(args.out_dir) / "similarity_matrix.csv"
            if not sim_path.exists():
                raise SystemExit("similarity_matrix.csv not found. Run --step compare first.")
            sim_df = pd.read_csv(sim_path, index_col=0)
        plot_similarity_heatmap(sim_df, args.out_dir, name="heatmap.png")

        # also drop clustermaps from summary if we have it
        if summary_df is not None:
            plot_clustermap(summary_df, args.out_dir, proportions=False)
            plot_clustermap(summary_df, args.out_dir, proportions=True)

    # SINGLE-FILE utilities
    if args.step == 'single-file':
        if not args.vector_file:
            raise SystemExit("--vector-file is required for step 'single-file'")
        single_file_plots(args.vector_file, args.out_dir)
        # optional: compare the single-file vectors to comparison vectors
        if args.comparison_dir or args.comparison_csv:
            single_file_compare(args.vector_file, args.comparison_dir or "", args.out_dir,
                                comparison_csv=args.comparison_csv)


if __name__ == '__main__':
    main()
