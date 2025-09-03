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

# Use 'X' consistently for stop in column names
AA_1 = list("ARNDCQEGHILKMFPSTWYV") + ["X"]                 # 21 symbols with X=stop
PAIR_COLS = [a + b for a in AA_1 for b in AA_1]             # 441 AA-pair columns, canonical order

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

def _to_one_letter_or_X(x: object) -> Optional[str]:
    """Map labels to 1-letter AA; STOP/*/TER/Ter/X => 'X'."""
    if x is None:
        return None
    s = str(x).strip()
    if not s:
        return None
    su = s.upper()
    if su in ("*", "STOP", "TER", "TERM", "X"):
        return "X"
    if len(su) == 1 and su in AA_1:
        return su
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

CANON_COLS = PAIR_COLS
CANON_SET  = set(CANON_COLS)

# Ensure all HEADING codes are valid AA pairs from AA (incl. X)
_invalid = [c for c in CANON_COLS if len(c)!=2 or c[0] not in AA_1 or c[1] not in AA_1]
if _invalid:
    raise ValueError(f"HEADING has invalid codes: {_invalid[:5]}{'...' if len(_invalid)>5 else ''}")
   
# ---------- IO & Vectorization ----------

def _read_one_sample_csv(fp: Path) -> pd.Series:
    """
    Accept one of:
      (a) tall format: columns include ST/END/(Count|count)
      (b) 21x21 matrix: rows=AA, cols=AA (any orientation; first col may be row labels)
      (c) vector-per-sample: a single row whose columns are AA-pair codes (AE, AG, …)
          with optional ID and/or SUM columns.
    Returns a Series indexed by CANON_COLS with float counts.
    """

    try:
        df = pd.read_csv(fp, sep=None, engine="python", compression="infer")
    except Exception:
        df = pd.read_csv(fp, compression="infer")

    # (c) single-row vectors (wide)
    present_pairs = [c for c in CANON_COLS if c in df.columns]
    if len(present_pairs) >= max(10, int(0.05 * len(CANON_COLS))) and len(df) == 1:
        row = df.iloc[0]
        vec = pd.Series(0.0, index=CANON_COLS, dtype=float)
        for c in present_pairs:
            try: vec[c] = float(row[c])
            except Exception: pass
        return vec

    # (a) tall format
    lower = {c.lower(): c for c in df.columns}
    def _has(*names): return all(n in lower for n in names)
    if _has('st','end') or _has('fromaa','toaa') or _has('refaa','altaa') or _has('from','to') or _has('ref_aa','alt_aa'):
        for L,R in [('st','end'),('fromaa','toaa'),('from','to'),('refaa','altaa'),('ref_aa','alt_aa')]:
            if L in lower and R in lower:
                df = df.rename(columns={lower[L]:'ST', lower[R]:'END'})
                break
        cnt_col = lower.get('count') or 'Count'
        if cnt_col not in df.columns:
            df[cnt_col] = 1
        df['ST'] = df['ST'].apply(_to_one_letter_or_X)
        df['END'] = df['END'].apply(_to_one_letter_or_X)
        df = df[df['ST'].isin(AA_1) & df['END'].isin(AA_1)]
        g = df.groupby(['ST','END'])[cnt_col].sum()
        vec = pd.Series(0.0, index=CANON_COLS, dtype=float)
        for (s,e), v in g.items():
            code = f"{s}{e}"
            if code in CANON_COLS:
                vec[code] += float(v)
        return vec

    # (b) 21×21 matrix
    if df.shape[1] >= 2:
        df = df.set_index(df.columns[0])
    rows_norm = [_to_one_letter_or_X(x) for x in df.index]
    cols_norm = [_to_one_letter_or_X(x) for x in df.columns]
    row_valid = sum(1 for r in rows_norm if r in AA_1)
    col_valid = sum(1 for c in cols_norm if c in AA_1)

    mat = df.copy()
    if row_valid < col_valid:
        mat = mat.T
        rows_norm, cols_norm = cols_norm, rows_norm
    mat.index = rows_norm
    mat.columns = cols_norm
    mat = mat.loc[[r for r in mat.index if r in AA_1],
                  [c for c in mat.columns if c in AA_1]]
    if mat.empty:
        raise ValueError(f"{fp.name}: unrecognized format (not tall, vector, or 21×21).")
    mat = mat.apply(pd.to_numeric, errors='coerce').fillna(0)

    vec = pd.Series(0.0, index=CANON_COLS, dtype=float)
    for r in mat.index:
        for c in mat.columns:
            code = r + c
            if code in CANON_COLS:
                vec[code] = float(mat.at[r, c])
    return vec


def _read_combined_vectors_table(fp: Path) -> pd.DataFrame | None:
    """Read multi-row combined vectors table (columns are AA-pair codes; +optional ID/SUM)."""
    try:
        df = pd.read_csv(fp, sep=None, engine="python", compression="infer")
    except Exception:
        df = pd.read_csv(fp, compression="infer")
    present = [c for c in CANON_COLS if c in df.columns]
    if len(present) >= max(50, int(0.2 * len(CANON_COLS))) and len(df) >= 2:
        if 'ID' not in df.columns:
            df = df.rename(columns={df.columns[0]: 'ID'})
        out = df[['ID'] + present].copy()
        out[CANON_COLS] = out.reindex(columns=CANON_COLS, fill_value=0).apply(pd.to_numeric, errors='coerce').fillna(0.0)
        out['SUM'] = out[CANON_COLS].sum(axis=1)
        return out[['ID','SUM'] + CANON_COLS]
    return None


def summarize_dir(dir_path: str, manifest_path: str | None = None) -> pd.DataFrame:
    """
    Recursively ingest .csv/.CSV/.tsv/.TSV/.csv.gz/.tsv.gz files.
    Accept tall, matrices, single-row vectors, and combined tables.
    """
    dirp = Path(dir_path)
    patterns = ["*.csv","*.CSV","*.tsv","*.TSV","*.csv.gz","*.CSV.GZ","*.tsv.gz","*.TSV.GZ"]
    files = sorted({p for pat in patterns for p in dirp.rglob(pat)})
    logging.info(f"[summarize_dir] Scanning {dir_path} — {len(files)} candidate file(s)")

    rows = []
    for fp in files:
        # combined table?
        comb = _read_combined_vectors_table(fp)
        if comb is not None:
            logging.info(f"[summarize_dir] Combined table: {fp.name} (+{len(comb)} rows)")
            rows.append(comb)
            continue
        # per-sample file
        try:
            vec = _read_one_sample_csv(fp)
            sid = fp.stem
            one = pd.DataFrame([{'ID': sid, 'SUM': float(np.nansum(vec))} | {c: float(vec[c]) for c in CANON_COLS}])
            rows.append(one)
        except Exception as e:
            logging.debug(f"[summarize_dir] Skipping {fp}: {e}")

    if not rows:
        raise FileNotFoundError(f"No usable CSV/TSV files found in {dir_path}")

    df = pd.concat(rows, ignore_index=True)
    if 'ID' not in df.columns:
        df = df.rename(columns={df.columns[0]: 'ID'})
    df = df[['ID','SUM'] + CANON_COLS]
    if manifest_path:
        ids = pd.read_table(manifest_path, header=None).iloc[:,0].astype(str)
        df = df[df['ID'].astype(str).isin(set(ids))]
    return df

def summarize_observed(observed_dir: str, out_dir: str, manifest_path: Optional[str]) -> pd.DataFrame:
    df = summarize_dir(observed_dir, manifest_path)
    out_path = Path(out_dir)/"summary.csv"
    df.to_csv(out_path, index=False)
    logging.info(f"Saved summary -> {out_path}")
    return df

def load_comparison(comparison_dir: str,
                    comparison_csv: Optional[str] = None,
                    manifest_path: Optional[str] = None,
                    out_dir: Optional[str] = None) -> pd.DataFrame:
    """
    Load comparison vectors. If a combined CSV isn’t provided/found, build it from the directory.
    Returns a DataFrame indexed by ID (no SUM), aligned to CANON_COLS.
    """
    # 1) Explicit file wins
    if comparison_csv:
        path = Path(comparison_csv)
        if not path.exists():
            raise FileNotFoundError(f"--comparison-csv not found: {path}")
        logging.info(f"Using comparison file -> {path}")
        df = pd.read_csv(path)

    else:
        # 2) Try common filenames inside the folder
        if not comparison_dir:
            raise SystemExit("Provide --comparison-dir or --comparison-csv for comparison steps.")
        comp_dir = Path(comparison_dir)
        candidates = ("comparison_vectors.csv", "comparison_summary.csv", "summary.csv")
        path = next((comp_dir / name for name in candidates if (comp_dir / name).exists()), None)

        if path is not None:
            logging.info(f"Using comparison file -> {path}")
            df = pd.read_csv(path)
        else:
            # 3) Fall back: summarize the comparison directory itself
            logging.info(f"No combined comparison file found in {comp_dir}. Summarizing that directory …")
            df = summarize_dir(str(comp_dir), manifest_path)
            if out_dir:
                outp = Path(out_dir) / "comparison_summary.csv"
                df.to_csv(outp, index=False)
                logging.info(f"Wrote comparison summary -> {outp}")

    # Normalize & align
    if 'ID' not in df.columns:
        df = df.rename(columns={df.columns[0]: 'ID'})
    df = df.set_index('ID')
    df = df.drop(columns=['SUM'], errors='ignore')

    missing = [c for c in CANON_COLS if c not in df.columns]
    if missing:
        logging.warning(f"Comparison vectors missing {len(missing)} AA codes (filled with 0). "
                        f"Examples: {missing[:8]}")

    return df.reindex(columns=CANON_COLS, fill_value=0.0)


# ---------- Comparisons & Plots ----------

def compare_vectors(observed_summary: pd.DataFrame,
                    comparison_df: pd.DataFrame,
                    out_dir: str,
                    outfile: str = "similarity_matrix.csv") -> pd.DataFrame:
    obs = observed_summary.set_index('ID').drop(columns=['SUM'], errors='ignore')
    common = [c for c in obs.columns if c in comparison_df.columns]
    if not common:
        raise ValueError("No overlapping AA-pair columns.")
    X = obs[common].to_numpy(dtype=float)
    Y = comparison_df[common].to_numpy(dtype=float)
    Xn = X / (np.linalg.norm(X, axis=1, keepdims=True) + 1e-12)
    Yn = Y / (np.linalg.norm(Y, axis=1, keepdims=True) + 1e-12)
    sims = Xn @ Yn.T
    sim_df = pd.DataFrame(sims, index=obs.index, columns=comparison_df.index)
    out_path = Path(out_dir)/outfile
    sim_df.to_csv(out_path)
    logging.info(f"Saved similarity matrix -> {out_path}")
    return sim_df

def plot_similarity_heatmap(sim_df: pd.DataFrame, out_dir: str, name: str = "heatmap.png"):
    plt.figure(figsize=(10,10))
    sb.heatmap(sim_df.astype(float), cmap='viridis')
    out_path = Path(out_dir)/name
    plt.tight_layout(); plt.savefig(out_path); plt.close()
    logging.info(f"Saved heatmap -> {out_path}")

def plot_clustermap_vectors(vectors_df, out_dir, normalize_rows=True, name=None, metric='cosine'):
    """
    Cluster map of sample x AA-pair vectors (optionally row-normalized to proportions).
    Drops rows/cols that would create NaN/Inf distances (e.g., all-zero rows for cosine).
    """
    X = vectors_df.set_index('ID').drop(columns=['SUM'], errors='ignore')
    # force numeric & clean non-finites
    X = X.apply(pd.to_numeric, errors='coerce').replace([np.inf, -np.inf], np.nan).fillna(0.0)

    # drop all-zero rows before cosine (undefined)
    if normalize_rows:
        row_sums = X.sum(axis=1)
        zero_rows = (row_sums == 0)
        if zero_rows.any():
            logging.warning(f"Clustermap: dropping {int(zero_rows.sum())} all-zero rows before normalization.")
            X = X.loc[~zero_rows]
        # proportions in %
        X = X.div(X.sum(axis=1), axis=0).fillna(0.0) * 100.0
    else:
        if metric == 'cosine':
            norms = np.linalg.norm(X.values, axis=1)
            zero_norm = (norms == 0)
            if zero_norm.any():
                logging.warning(f"Clustermap: dropping {int(zero_norm.sum())} rows with zero norm for cosine metric.")
                X = X.loc[~zero_norm]

    # drop all-zero columns (no information)
    zero_cols = (X.sum(axis=0) == 0)
    if zero_cols.any():
        logging.info(f"Clustermap: dropping {int(zero_cols.sum())} all-zero columns.")
        X = X.loc[:, ~zero_cols]

    if X.shape[0] < 2 or X.shape[1] < 2:
        logging.warning("Clustermap: not enough data to cluster after filtering; skipping.")
        return

    if name is None:
        name = "clustermap_proportions.png" if normalize_rows else "clustermap_counts.png"
    out_path = Path(out_dir) / name

    try:
        g = sb.clustermap(X, metric=metric, cmap='viridis', figsize=(14, 12))
    except ValueError as e:
        logging.warning(f"Clustermap failed with metric={metric} ({e}); retrying with 'euclidean'.")
        g = sb.clustermap(X, metric='euclidean', cmap='viridis', figsize=(14, 12))

    g.fig.savefig(out_path)
    plt.close(g.fig)
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
    # keep only known pair columns
    cols_present = [c for c in CANON_COLS if c in df.columns]
    if 'SUM' not in df.columns:
        df['SUM'] = df[cols_present].sum(axis=1)
    return df[['ID','SUM'] + cols_present]

def plot_single_file(vector_csv: str, out_dir: str):
    df = _read_summary_csv(vector_csv)
    df = df.set_index('ID')
    vec = df.drop(columns=['SUM'], errors='ignore').iloc[0]
    # Counts
    plt.figure(figsize=(18, 3))
    sb.heatmap([vec.values], cmap='viridis', cbar_kws={'label':'Count'})
    plt.yticks([0], [df.index[0]])
    plt.xticks(ticks=np.arange(len(vec)), labels=vec.index, rotation=90)
    out_counts = Path(out_dir)/"single_file_counts.png"
    plt.tight_layout(); plt.savefig(out_counts); plt.close()
    # Proportions
    prop = (vec / (vec.sum() or 1)) * 100.0
    plt.figure(figsize=(18, 3))
    sb.heatmap([prop.values], cmap='viridis', cbar_kws={'label':'% of substitutions'})
    plt.yticks([0], [df.index[0]])
    plt.xticks(ticks=np.arange(len(prop)), labels=prop.index, rotation=90)
    out_props = Path(out_dir)/"single_file_proportions.png"
    plt.tight_layout(); plt.savefig(out_props); plt.close()
    logging.info(f"Saved single-file plots -> {out_counts}, {out_props}")

def single_file_compare(vector_file: str, comparison_dir: str, out_dir: str, comparison_csv: Optional[str] = None):
    """
    Compare vectors from a precomputed summary file to comparison vectors.
    Saves single_similarity.csv and single_similarity_heatmap.png
    """
    df = _read_summary_csv(vector_file)
    comp = load_comparison(comparison_dir, comparison_csv=comparison_csv)
    sim = compare_vectors(df, comp, out_dir, outfile="single_similarity.csv")
    plot_similarity_heatmap(sim, out_dir, name="single_similarity_heatmap.png")

def compare_matrices(mat_a_path, mat_b_path, out_dir):
    """
    Compare two 21×21 AA matrices:
      • align and save counts for A and B
      • compute A−B proportion difference and save plain + clustered heatmaps
      • compute per-cell Fisher tests with BH-FDR
    Outputs are written under <out_dir>/matrices/.
    """
    from pathlib import Path
    import logging
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt

    # Optional deps
    try:
        from scipy.stats import fisher_exact
    except Exception:
        fisher_exact = None
    try:
        from scipy.cluster.hierarchy import linkage, leaves_list
    except Exception:
        linkage = None
        leaves_list = None

    mats = Path(out_dir) / "matrices"
    mats.mkdir(parents=True, exist_ok=True)

    # Load and align
    A = pd.read_csv(mat_a_path, index_col=0)
    B = pd.read_csv(mat_b_path, index_col=0)
    idx = sorted(set(A.index) | set(B.index))
    cols = sorted(set(A.columns) | set(B.columns))
    A = A.reindex(index=idx, columns=cols, fill_value=0)
    B = B.reindex(index=idx, columns=cols, fill_value=0)

    # Save aligned counts
    A.to_csv(mats / "A_counts.csv")
    B.to_csv(mats / "B_counts.csv")

    # Proportions & enrichment (A − B)
    At = max(int(A.values.sum()), 1)
    Bt = max(int(B.values.sum()), 1)
    Ap = A / At
    Bp = B / Bt
    D = Ap - Bp
    D.to_csv(mats / "A_minus_B_enrichment.csv")

    # Plain (unclustered) heatmap of A−B
    vmax = float(np.nanmax(np.abs(D.values))) if np.isfinite(np.nanmax(np.abs(D.values))) else 1.0
    if vmax <= 0:
        vmax = 1.0  # avoid zero range
    fig, ax = plt.subplots(figsize=(8, 6))
    im = ax.imshow(D.values, aspect="auto", vmin=-vmax, vmax=vmax)
    fig.colorbar(im, ax=ax)
    ax.set_xticks(range(len(D.columns)))
    ax.set_xticklabels(D.columns, rotation=90, fontsize=8)
    ax.set_yticks(range(len(D.index)))
    ax.set_yticklabels(D.index, fontsize=8)
    ax.set_title("A − B enrichment (proportion)")
    fig.tight_layout()
    fig.savefig(mats / "A_minus_B_enrichment_heatmap.png", dpi=150)
    plt.close(fig)

    # Clustered heatmap of A−B (Ward; rows + cols)
    if linkage is not None and leaves_list is not None:
        try:
            row_Z = linkage(D.values, method="ward", optimal_ordering=True)
            col_Z = linkage(D.values.T, method="ward", optimal_ordering=True)
            r_ord = leaves_list(row_Z)
            c_ord = leaves_list(col_Z)
            D_ord = D.iloc[list(r_ord), list(c_ord)]

            # Save clustered CSV + orders
            D_ord.to_csv(mats / "A_minus_B_enrichment.clustered.csv")
            pd.Series(D_ord.index, name="row_order").to_csv(mats / "A_minus_B_row_order.csv", index=False)
            pd.Series(D_ord.columns, name="col_order").to_csv(mats / "A_minus_B_col_order.csv", index=False)

            # Plot clustered heatmap
            fig2, ax2 = plt.subplots(figsize=(8, 6))
            im2 = ax2.imshow(D_ord.values, aspect="auto", vmin=-vmax, vmax=vmax)
            fig2.colorbar(im2, ax=ax2)
            ax2.set_xticks(range(len(D_ord.columns)))
            ax2.set_xticklabels(D_ord.columns, rotation=90, fontsize=8)
            ax2.set_yticks(range(len(D_ord.index)))
            ax2.set_yticklabels(D_ord.index, fontsize=8)
            ax2.set_title("A − B enrichment (clustered)")
            fig2.tight_layout()
            fig2.savefig(mats / "A_minus_B_enrichment_clustered.png", dpi=150)
            plt.close(fig2)
        except Exception as e:
            logging.warning(f"[matrix-compare] clustering failed: {e}")

    # Per-cell stats (Fisher + BH-FDR)
    rows = []
    for i in idx:
        for j in cols:
            a = int(A.loc[i, j])        # cell in A
            b = int(B.loc[i, j])        # cell in B
            c = int(At - a)             # all other cells in A
            d = int(Bt - b)             # all other cells in B

            # Haldane–Anscombe correction for zeros
            a2, b2, c2, d2 = a + 0.5, b + 0.5, c + 0.5, d + 0.5
            if fisher_exact is not None:
                try:
                    _, p = fisher_exact([[a2, c2], [b2, d2]], alternative="two-sided")
                except Exception:
                    p = 1.0
            else:
                p = 1.0

            rows.append({
                "AA_from": i,
                "AA_to": j,
                "A_count": a,
                "B_count": b,
                "A_prop": a / At,
                "B_prop": b / Bt,
                "log2_enrichment": float(np.log2((a / At + 1e-12) / (b / Bt + 1e-12))),
                "pval": float(p),
            })

    df = pd.DataFrame(rows)
    pvals = df["pval"].values
    n = len(pvals)
    if n:
        order = np.argsort(pvals)
        ranks = np.arange(1, n + 1)
        q = np.empty_like(pvals, dtype=float)
        q[order] = np.minimum.accumulate((pvals[order] * n / ranks)[::-1])[::-1]
        df["qval"] = np.clip(q, 0.0, 1.0)
    else:
        df["qval"] = []

    df.to_csv(mats / "cell_stats.csv", index=False)

    return {
        "A_counts": str((mats / "A_counts.csv").resolve()),
        "B_counts": str((mats / "B_counts.csv").resolve()),
        "A_minus_B_enrichment": str((mats / "A_minus_B_enrichment.csv").resolve()),
        "cell_stats": str((mats / "cell_stats.csv").resolve()),
    }

# ---------- CLI ----------

def parse_args():
    p = argparse.ArgumentParser(description="Comparison Pipeline")
    p.add_argument('--observed-dir', '--observed_dir', '--obs-dir', '--obs_dir',
                   dest='observed_dir',
                   help='Directory of observed AA CSVs (tall files or 21x21 matrices)')
    p.add_argument('--comparison-dir', '--comparison_dir', '--comp-dir', '--comp_dir',
                   dest='comparison_dir',
                   help='Directory containing comparison vectors (or per-sample CSVs to be summarized)')
    p.add_argument('--comparison-csv', '--comparison_csv', dest='comparison_csv',
                   help='Optional explicit path to comparison/reference combined CSV')
    p.add_argument('--manifest', help='Optional manifest (first column are sample IDs) to filter observed')
    p.add_argument('--out-dir', '--out_dir', dest='out_dir', required=True,
                   help='Directory for outputs')
    p.add_argument('--vector-file', '--vector_file', dest='vector_file',
                   help='Path to an existing summary CSV (to skip summarize or for single-file)')
    p.add_argument('--step', choices=['summarize','compare','heatmap','single-file', 'matrix-compare','all'],
                   default='all', help='Which step(s) to run')
    p.add_argument('--cluster-metric', choices=['cosine','euclidean','correlation'],
                   default='cosine', help='Distance metric for clustermap (default: cosine)')
    p.add_argument('--matrix-a', help='21x21 AA matrix CSV for set A')
    p.add_argument('--matrix-b', help='21x21 AA matrix CSV for set B')

    return p.parse_args()

# ---------- Main ----------

def main():
    args = parse_args()
    out_dir = Path(args.out_dir); out_dir.mkdir(parents=True, exist_ok=True)

    # EARLY EXIT: matrix-compare does not need comparison vectors
    if args.step == 'matrix-compare':
        if not args.matrix_a or not args.matrix_b:
            raise SystemExit("--matrix-a and --matrix-b are required for 'matrix-compare'")
        out_paths = compare_matrices(args.matrix_a, args.matrix_b, args.out_dir)
        logging.info(f"[matrix-compare] wrote: {out_paths}")
        return 0

    summary_df = None
    sim_df = None

    # SUMMARIZE (or load existing)
    if args.step in ('summarize', 'all'):
        if not args.observed_dir:
            raise SystemExit("--observed-dir is required for step 'summarize' or 'all'")
        summary_df = summarize_observed(args.observed_dir, args.out_dir, args.manifest)
        # clustermaps for observed cohort
        plot_clustermap_vectors(summary_df.copy(), args.out_dir,
                                normalize_rows=False, name="clustermap_counts.png",
                                metric=args.cluster_metric)
        plot_clustermap_vectors(summary_df.copy(), args.out_dir,
                                normalize_rows=True, name="clustermap_proportions.png",
                                metric=args.cluster_metric)
    else:
        # Try to get a summary from vector_file or out_dir/summary.csv
        if args.vector_file and Path(args.vector_file).exists():
            summary_df = _read_summary_csv(args.vector_file)
        else:
            vf = Path(args.out_dir) / 'summary.csv'
            if vf.exists():
                summary_df = _read_summary_csv(str(vf))

    # COMPARE
    if args.step in ('compare', 'all'):
        if summary_df is None:
            raise SystemExit("No summary available. Run with --step summarize or provide --vector-file/summary.csv.")
        if not args.comparison_dir and not args.comparison_csv:
            raise SystemExit("--comparison-dir or --comparison-csv is required for 'compare'/'all'")
        comp = load_comparison(args.comparison_dir or "",
                               comparison_csv=args.comparison_csv,
                               manifest_path=args.manifest,
                               out_dir=args.out_dir)
        sim_df = compare_vectors(summary_df, comp, args.out_dir)
        plot_similarity_heatmap(sim_df, args.out_dir, name="heatmap.png")

    # HEATMAP only (re-use existing similarity)
    if args.step == 'heatmap':
        sim_path = Path(args.out_dir) / "similarity_matrix.csv"
        if not sim_path.exists():
            raise SystemExit("similarity_matrix.csv not found. Run --step compare first.")
        sim_df = pd.read_csv(sim_path, index_col=0)
        plot_similarity_heatmap(sim_df, args.out_dir, name="heatmap.png")
        if summary_df is not None:
            plot_clustermap_vectors(summary_df.copy(), args.out_dir, normalize_rows=False, name="clustermap_counts.png")
            plot_clustermap_vectors(summary_df.copy(), args.out_dir, normalize_rows=True, name="clustermap_proportions.png")

    # SINGLE-FILE utilities
    if args.step == 'single-file':
        if not args.vector_file:
            raise SystemExit("--vector-file is required for step 'single-file'")
        plot_single_file(args.vector_file, args.out_dir)
        if args.comparison_dir or args.comparison_csv:
            single_file_compare(args.vector_file, args.comparison_dir or "", args.out_dir,
                                comparison_csv=args.comparison_csv)


   

if __name__ == '__main__':
    main()
   
