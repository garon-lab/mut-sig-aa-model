#!/usr/bin/env python3
from __future__ import annotations
"""
Signature-Based AA Variant Modeling (Vectorized)

This script takes a 12-element vector or COSMIC mutational signature (context proportions) and produces
expected amino acid substitution profiles based on embedded base matrices.

Features:
1. Input CSV: ID + 12 signature proportions (AC, AG, AT, CA, CG, CT, GA, GC, GT, TA, TC, TG) OR COSMIC signature table with flag.
2. Base substitution matrices are embedded as NumPy arrays—no external files needed.
3. Computes weighted sums, scaling, and row-wise normalization to amino-acid frequency targets using vectorized operations.
4. Flattens 21×21 expected matrices into vectors for downstream analysis.
5. Optionally plots the expected matrix as a heatmap with customizable appearance.

Usage:
    python signature_modeling.py \
        --signature_vector <CSV of signature proportions> \
        --out_dir <output directory> \
        [--step model|heatmap|both] \
        [--log_level DEBUG|INFO|WARNING|ERROR]

Arguments:
    --signature_vector   CSV with columns: ID, AC, AG, AT, CA, CG, CT, GA, GC, GT, TA, TC, TG or COSMIC signature table.
    --out_dir            Directory to write expected vectors and heatmap
    --step               Step to run: model, heatmap, or both (default: both)
    --log_level          Logging level (default: INFO)

Dependencies:
    pandas, numpy, matplotlib, seaborn
"""

import argparse
import logging
from pathlib import Path

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import re
import seaborn as sb

try:
    import scipy
except ImportError:
    logging.error("SciPy is required for clustering; install with `pip install scipy`.")

_SIX = ["C>A","C>G","C>T","T>A","T>C","T>G"]
_TWELVE = ["AC","AG","AT","CA","CG","CT","GA","GC","GT","TA","TC","TG"]
# Map 6 to 12 by distributing the pyrimidine-centered classes:
_SIX_TO_TWELVE = {
    "C>A": ["CA","GT"],
    "C>G": ["CG","GC"],
    "C>T": ["CT","GA"],
    "T>A": ["TA","AT"],
    "T>C": ["TC","AG"],
    "T>G": ["TG","AC"],
}

def cosmic96_to_12(s96: pd.Series) -> pd.Series:
    """
    s96: Series indexed by 96 'Type' contexts like A[C>A]A with values.
    Returns a 12-channel Series (AC..TG) normalized to sum==1 if input sums to 1.
    """
    # collapse 96 -> 6 by extracting the middle base change inside [X>Y]
    six = pd.Series(0.0, index=_SIX, dtype=float)
    for ctx, v in s96.items():
        m = re.search(r"\[([ACGT]>\w)\]", ctx)
        if not m:
            continue
        mid = m.group(1).upper()        # e.g., C>A
        if mid in six.index:
            six[mid] += float(v) if pd.notna(v) else 0.0
    # distribute 6 -> 12 by splitting evenly into the paired 12 slots
    twelve = pd.Series(0.0, index=_TWELVE, dtype=float)
    for k, dests in _SIX_TO_TWELVE.items():
        portion = six[k] / len(dests) if len(dests) else 0.0
        for d in dests:
            twelve[d] += portion
    return twelve
 
# Plot customization constants
def _parse_figsize(s: str):
    s = (s or "12x8").lower().replace(" ", "").replace(",", "x")
    try:
        w, h = s.split("x")
        return (float(w), float(h))
    except Exception:
        logging.warning(f"Bad --figsize '{s}', falling back to 12x8")
        return (12.0, 8.0)

FIGSIZE = (12, 8)
COLORMAP = 'viridis'
HEATMAP_KWARGS = {
    'cmap': COLORMAP,
    'annot': False,
    'linewidths': 0.5,
    'linecolor': 'gray'
}

# Logging initialization
def init_logging(level_str):
    level = getattr(logging, level_str.upper(), None)
    if not isinstance(level, int):
        raise ValueError(f"Invalid log level: {level_str}")
    logging.basicConfig(level=level, format='%(asctime)s %(levelname)s: %(message)s')

# Helpers
def validate_context_df(df):
    cols = set(df.columns)
    extras = [c for c in cols - set(CONTEXTS) if c.lower() not in {"id","sample","sample_id"}]
    missing = [c for c in CONTEXTS if c not in cols]
    if missing or extras:
        msg = []
        if missing: msg.append(f"Missing: {', '.join(missing)}")
        if extras:  msg.append(f"Unexpected: {', '.join(extras)}")
        raise ValueError("Signature validation error: " + "; ".join(msg))
    logging.info("Signature vector validated.")

def load_cosmic_signature_table(path: str | Path) -> tuple[str, pd.Series]:
    """
    Load a COSMIC-style signature .txt file where the first row has headers like:
        Type <TAB> SBS4
    and subsequent rows look like:
        A[C>A]A <TAB> 0.042...
    Returns (signature_name, Series(index=context, values=float)).
    """
    path = Path(path)
    df = pd.read_csv(path, sep=r"\t", engine="python", dtype=str, comment="#").rename(
        columns=lambda c: c.strip()
    )
    # Identify the single non-'Type' column
    cols = [c for c in df.columns if c.lower() != "type"]
    if len(cols) != 1:
        # Fall back to filename stem if odd formatting
        sig_name = path.stem
        # If there are multiple non-Type columns, just pick the last one
        sig_col = cols[-1] if cols else None
    else:
        sig_name = cols[0]
        sig_col  = cols[0]
    if sig_col is None or "Type" not in df.columns:
        raise ValueError(f"Could not parse signature file: {path}")
    df["Type"] = df["Type"].astype(str).str.strip()
    s = pd.to_numeric(df[sig_col], errors="coerce")
    sig = pd.Series(s.values, index=df["Type"].values, dtype=float).dropna()
    return sig_name, sig

def load_signatures_from_path(p: str | Path) -> dict[str, pd.Series]:
    """
    If 'p' is a file: load that one signature.
    If 'p' is a directory: load all *.txt files inside.
    Returns dict[name] -> Series(context->value).
    """
    p = Path(p)
    sigs: dict[str, pd.Series] = {}
    files = [p] if p.is_file() else sorted(p.glob("*.txt"))
    if not files:
        raise FileNotFoundError(f"No signature file(s) found at: {p}")
    for f in files:
        name, sig = load_cosmic_signature_table(f)
        # If duplicate names, disambiguate with filename
        if name in sigs:
            name = f.stem
        sigs[name] = sig
    return sigs

def _normalize_nonnegative(x: pd.Series | np.ndarray, eps: float = 1e-12) -> np.ndarray:
    arr = np.asarray(x, dtype=float)
    arr[np.isnan(arr)] = 0.0
    arr[arr < 0]      = 0.0
    s = arr.sum()
    return arr / (s + eps)

def _cosine_sim(a: np.ndarray, b: np.ndarray, eps: float = 1e-12) -> float:
    num = float(np.dot(a, b))
    den = float(np.linalg.norm(a) * np.linalg.norm(b)) + eps
    return num / den

def _pearson(a: np.ndarray, b: np.ndarray, eps: float = 1e-12) -> float:
    a = a - a.mean()
    b = b - b.mean()
    num = float(np.dot(a, b))
    den = float(np.linalg.norm(a) * np.linalg.norm(b)) + eps
    return num / den

def _l1(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.sum(np.abs(a - b)))

def _l2(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.linalg.norm(a - b))

def compare_signatures_heatmap(
    observed: pd.DataFrame,
    signatures: dict[str, pd.Series],
    out_png: str | Path,
    out_csv: str | Path,
    *,
    metric: str = "cosine",
    normalize_observed: bool = True,
    normalize_signatures: bool = True,
    title: str = "Observed vs Signatures"
) -> pd.DataFrame:
    """
    Compare observed vectors (rows = samples, columns = contexts like A[C>A]A)
    against signatures (Series, index=contexts). Saves heatmap PNG + CSV.
    Returns the score matrix DataFrame (n_samples x n_signatures).
    """
    # Build full context set
    all_contexts = set(observed.columns)
    for sig in signatures.values():
        all_contexts |= set(sig.index)
    contexts = sorted(all_contexts)

    obs = observed.copy().reindex(columns=contexts, fill_value=0.0)

    sig_mat = {}
    for name, sig in signatures.items():
        s = sig.reindex(index=contexts, fill_value=0.0)
        if normalize_signatures:
            s = pd.Series(_normalize_nonnegative(s), index=s.index)
        sig_mat[name] = s.values
    sig_names = list(sig_mat.keys())
    S = np.vstack([sig_mat[n] for n in sig_names])  # (k, d)

    X = obs.values.astype(float)
    if normalize_observed:
        row_sums = X.sum(axis=1, keepdims=True) + 1e-12
        X = np.clip(X, 0, None) / row_sums

    m = metric.lower()
    if m == "cosine":
        scorer, is_similarity = _cosine_sim, True
    elif m == "pearson":
        scorer, is_similarity = _pearson, True
    elif m == "l1":
        scorer, is_similarity = _l1, False
    elif m == "l2":
        scorer, is_similarity = _l2, False
    else:
        raise ValueError("metric must be one of: cosine, pearson, l1, l2")

    scores = np.zeros((X.shape[0], S.shape[0]), dtype=float)
    for i in range(X.shape[0]):
        xi = X[i]
        for j in range(S.shape[0]):
            sj = S[j]
            scores[i, j] = scorer(xi, sj)

    out_df = pd.DataFrame(scores, index=obs.index, columns=sig_names)
    out_csv = Path(out_csv); out_df.to_csv(out_csv)

    # Heatmap (closed properly to avoid >20 open figs)
    fig, ax = plt.subplots(figsize=(max(6, 0.5*len(sig_names)+2), max(4, 0.4*len(obs.index)+2)))
    try:
        im = ax.imshow(out_df.values, aspect="auto", interpolation="nearest",
                       cmap=("viridis" if is_similarity else "magma"))
        ax.set_xticks(range(len(sig_names))); ax.set_xticklabels(sig_names, rotation=45, ha="right")
        ax.set_yticks(range(len(out_df.index))); ax.set_yticklabels(out_df.index)
        ax.set_xlabel("Signatures"); ax.set_ylabel("Observed (samples)")
        ax.set_title(f"{title} — {metric}")
        cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label("similarity" if is_similarity else "distance")
        fig.tight_layout()
        out_png = Path(out_png); fig.savefig(out_png, dpi=150, bbox_inches="tight")
    finally:
        plt.close(fig)

    return out_df
 
# ---- helpers for plotting CLI ----
AA_ORDER = ["A","R","N","D","C","Q","E","G","H","I",
            "L","K","M","F","P","S","T","W","Y","V","*"]

def _validate_and_coerce(name, M):
    # Ensure list-of-lists and coerce each cell to float
    if not isinstance(M, (list, tuple)):
        raise ValueError(f"{name} must be a list of rows, got {type(M)}")

    # Strip strings, handle ints/floats, and check row lengths
    coerced = []
    row_lengths = []
    for r_idx, row in enumerate(M):
        if not isinstance(row, (list, tuple)):
            raise ValueError(f"{name} row {r_idx} is not a list/tuple (type {type(row)})")
        new_row = []
        for c_idx, val in enumerate(row):
            try:
                # common issue: embedded strings like ' 0 ' or '0.0'
                new_row.append(float(str(val).strip()))
            except Exception as e:
                raise ValueError(
                    f"{name} has a non-numeric cell at row {r_idx}, col {c_idx}: {val!r}"
                ) from e
        row_lengths.append(len(new_row))
        coerced.append(new_row)

    # Check 21 rows
    if len(coerced) != 21:
        raise ValueError(f"{name} must have 21 rows (one per AA in {AA_ORDER}), found {len(coerced)}")

    # Check every row has 21 columns
    bad = [(i, L) for i, L in enumerate(row_lengths) if L != 21]
    if bad:
        sample = ", ".join([f"row{i}={L}" for i, L in bad[:5]])
        raise ValueError(
            f"{name} must be 21×21; mismatched row lengths found ({sample}{' ...' if len(bad)>5 else ''})."
        )

    import numpy as np
    return np.array(coerced, dtype=float)
 
# Amino-acid row targets
ROW_TARGETS = np.array([
    1.444267, 0.459934, 1.009417, 1.503128, 0.754421,
    1.352014, 0.542943, 0.915879, 1.217245, 2.061498,
    0.464837, 0.758762, 1.298454, 0.999901, 1.181556,
    1.741092, 1.117519, 1.253110, 0.317910, 0.554111,
    0.055119  # STOP
])

# Mutation contexts
CONTEXTS = ['AC','AG','AT','CA','CG','CT','GA','GC','GT','TA','TC','TG']
TOTAL_DIVISOR = 438.0
SCALE_FACTOR = 13.0  # scale = 1 + number of contexts

# Amino acid labels and SUB_VECTOR labels
# Replace your AA_LABELS / SUB_VECTOR block with this:
AA_LABELS = ["A","R","N","D","C","Q","E","G","H","I",
             "L","K","M","F","P","S","T","W","Y","V","*"]
SUB_VECTOR = ['ID'] + [f"{r}{c}" for r in AA_LABELS for c in AA_LABELS]


# Embedded base 21x21 matrices:
AC = [
 [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[2,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,3,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,1,0,2,0,0,2,0,0,0,0],
[0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,2,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,1,0,0,0,0,0,2,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,4,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0],
[0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,1,0]
]

AG = [
 [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,3,0,0,0],
[0,0,0,2,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0],
[0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0],
[0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0]
]

AT = [
 [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0],
[0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,2,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,1,0,0,1,1,0,0,0,0,0,0,0,0,2],
[0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,2,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,1,0,0,2,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,1],
[0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,1,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,1,0]
]

CA = [
 [0,0,2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1],
[0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,2,0,1,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,3,0,0,1,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,2,0,0,0,0,0,0,2,0,0,4,0,0,0,0],
[0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,2,2],
[0,0,0,0,0,0,0,0,2,0,0,2,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
]

CG = [
 [0,0,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0],
[0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,2,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0],
[4,0,0,0,0,0,0,0,0,0,0,0,0,0,4,0,0,0,0,0,0],
[0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,2,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,1],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,2,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
]

CT = [
 [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,4,0,0,0,0,0,4,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2],
[0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,2],
[0,0,0,0,2,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,3,0,0,1,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
]

GA = [
 [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0],
[0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,2,2,0,0,0,0,0,0,0,0,0,0,2,2,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,2,0,2,0,0,0,0,2,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,3,0,0,1,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
]

GC = [
 [0,0,0,0,0,0,0,0,0,0,0,0,4,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0],
[0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,1,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[4,0,0,0,0,0,0,0,0,0,0,0,0,0,4,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0],
[0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,4,0,0,1,2,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0],
[0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0]
]

GT = [
 [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4,0,0,0,0,0],
[0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0],
[0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4,1,0,1],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0],
[0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,1,0,4,1,0,0,0,0,1,0,0,0,0,0],
[0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,2,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0],
[0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0]
]

TA = [
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,2],
[0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,2,0,1,0,0,0,0,0,0,0,0,0,2,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,1,0,0,2,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,2,1,0,0,1,0,0,2,0,0,0,0,0,0,2],
[0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,4,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,2],
[0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,1],
[0,0,0,0,0,0,0,0,2,0,0,0,0,0,2,0,0,0,0,0,0]
]

TC = [
 [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,2,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,4,0,0,2,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,4,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0],
[0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,2,1,0,0,0,0,0,0]
]

TG = [
 [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0],
[0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,2,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,2,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,2,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,4,2,0,2,1,0,1],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1],
[0,0,0,2,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
]


base_matrices = {
    "AC": _validate_and_coerce("AC", AC),
    "AG": _validate_and_coerce("AG", AG),
    "AT": _validate_and_coerce("AT", AT),
    "CA": _validate_and_coerce("CA", CA),
    "CG": _validate_and_coerce("CG", CG),
    "CT": _validate_and_coerce("CT", CT),
    "GA": _validate_and_coerce("GA", GA),
    "GC": _validate_and_coerce("GC", GC),
    "GT": _validate_and_coerce("GT", GT),
    "TA": _validate_and_coerce("TA", TA),
    "TC": _validate_and_coerce("TC", TC),
    "TG": _validate_and_coerce("TG", TG),
}

# ---- Context order for the 12 signature channels (must match your CSV columns) ----
CONTEXT_ORDER = ["AC","AG","AT","CA","CG","CT","GA","GC","GT","TA","TC","TG"]

# Sanity: make sure all 12 are present
missing = [k for k in CONTEXT_ORDER if k not in base_matrices]
if missing:
    raise KeyError(f"Missing matrices for: {missing}")

# Each base matrix is 21x21; stack them into a (12, 21, 21) tensor
BASE_STACK = np.stack([base_matrices[k] for k in CONTEXT_ORDER], axis=0)
if BASE_STACK.shape != (12, 21, 21):
    raise ValueError(f"BASE_STACK shape is {BASE_STACK.shape}, expected (12, 21, 21)")

def build_expected(df):
    """Vectorized expected matrix computation."""
    # weights shape: (N, 12) using the same CONTEXT_ORDER as BASE_STACK
    weights = df[CONTEXT_ORDER].to_numpy(dtype=float)  # rows: samples, cols: AC..TG
    # Expected tensors: (N, 21, 21)
    mats = np.tensordot(weights, BASE_STACK, axes=([1], [0]))

    mats /= TOTAL_DIVISOR
    row_sums = mats.sum(axis=2)
    norm_factors = (ROW_TARGETS * len(AA_LABELS)) / np.where(row_sums==0, 1, row_sums)
    mats = mats * norm_factors[:,:,None]
    flat = mats.reshape(mats.shape[0], -1)
    return pd.DataFrame(flat, index=df.index, columns=SUB_VECTOR[1:])


def plot_heatmap(exp_df, out_dir, *, cmap, figsize, annot, linewidths, linecolor, vmin, vmax):
    """Plot heatmap of expected substitution profiles (no clustering)."""
    fig, ax = plt.subplots(figsize=figsize)
    sb.heatmap(
        exp_df.astype(float),
        ax=ax,
        cmap=cmap,
        annot=annot,
        linewidths=linewidths,
        linecolor=linecolor,
        vmin=vmin,
        vmax=vmax,
    )
    plt.tight_layout()
    path = Path(out_dir) / 'expected_heatmap.png'
    try:
        fig.savefig(path, dpi=150, bbox_inches="tight")
    finally:
        plt.close(fig)
    logging.info(f"Saved heatmap to {path}")

def plot_heatmap_clustered(
    exp_df,
    out_dir,
    *,
    cmap="viridis",
    figsize=(12,8),
    vmin=None,
    vmax=None,
    cluster="both",                 # 'none','rows','cols','both'
    method="average",               # linkage
    metric="euclidean",             # distance metric
    standardize="none"              # 'none','row','col'
):
    """
    Clustered heatmap using seaborn.clustermap with guards:
    - If rows<2, row clustering is disabled.
    - If cols<2, column clustering is disabled.
    - If both disabled, fall back to a plain heatmap.
    Always writes:
      - expected_heatmap_clustered.png
      - cluster_row_order.csv
      - cluster_col_order.csv
    """
    out_dir = Path(out_dir)
    n_rows, n_cols = exp_df.shape
    req_row = cluster in ("rows","both")
    req_col = cluster in ("cols","both")
    row_cluster = req_row and (n_rows >= 2)
    col_cluster = req_col and (n_cols >= 2)

    if not row_cluster and req_row:
        logging.warning("[cluster] Not clustering rows (need >=2 rows).")
    if not col_cluster and req_col:
        logging.warning("[cluster] Not clustering cols (need >=2 cols).")

    # If nothing can be clustered, write a simple heatmap + identity orders
    if not row_cluster and not col_cluster and cluster != "none":
        logging.warning("[cluster] Neither axis can be clustered; falling back to non-clustered heatmap.")
        fig, ax = plt.subplots(figsize=figsize)
        try:
            sb.heatmap(exp_df.astype(float), ax=ax, cmap=cmap, vmin=vmin, vmax=vmax)
            ax.set_xlabel("Substitution (21×21 flattened)")
            ax.set_ylabel("Samples")
            fig.tight_layout()
            (out_dir / "expected_heatmap_clustered.png").parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(out_dir / "expected_heatmap_clustered.png", dpi=150, bbox_inches="tight")
        finally:
            plt.close(fig)
        # Identity orders
        pd.DataFrame({"row_index": np.arange(n_rows), "row_label": exp_df.index.to_numpy()}) \
          .to_csv(out_dir / "cluster_row_order.csv", index=False)
        pd.DataFrame({"col_index": np.arange(n_cols), "col_label": exp_df.columns.to_numpy()}) \
          .to_csv(out_dir / "cluster_col_order.csv", index=False)
        logging.info("Saved cluster order to cluster_row_order.csv and cluster_col_order.csv")
        return

    # Map standardize -> seaborn params
    cg_kwargs = {}
    if standardize == "row":
        cg_kwargs["standard_scale"] = 1
    elif standardize == "col":
        cg_kwargs["standard_scale"] = 0

    # clustermap (ensure we close the fig to avoid figure leaks)
    g = sb.clustermap(
        exp_df.astype(float),
        cmap=cmap,
        vmin=vmin, vmax=vmax,
        method=method,
        metric=metric,
        row_cluster=row_cluster,
        col_cluster=col_cluster,
        figsize=figsize,
        **cg_kwargs
    )
    try:
        # Derive orders (fallback to identity if dendrogram missing)
        if row_cluster and hasattr(g, "dendrogram_row") and getattr(g.dendrogram_row, "reordered_ind", None) is not None:
            row_order = np.asarray(g.dendrogram_row.reordered_ind, dtype=int)
        else:
            row_order = np.arange(n_rows, dtype=int)

        if col_cluster and hasattr(g, "dendrogram_col") and getattr(g.dendrogram_col, "reordered_ind", None) is not None:
            col_order = np.asarray(g.dendrogram_col.reordered_ind, dtype=int)
        else:
            col_order = np.arange(n_cols, dtype=int)

        # Save cluster orders (labels mirrored to current df)
        pd.DataFrame({
            "row_index": row_order,
            "row_label": exp_df.index.to_numpy()[row_order]
        }).to_csv(out_dir / "cluster_row_order.csv", index=False)

        pd.DataFrame({
            "col_index": col_order,
            "col_label": exp_df.columns.to_numpy()[col_order]
        }).to_csv(out_dir / "cluster_col_order.csv", index=False)

        g.ax_heatmap.set_xlabel("Substitution (21×21 flattened)")
        g.ax_heatmap.set_ylabel("Samples")
        g.fig.tight_layout()
        g.fig.savefig(out_dir / "expected_heatmap_clustered.png", dpi=150, bbox_inches="tight")
        logging.info("Saved clustered heatmap to expected_heatmap_clustered.png")
        logging.info("Saved cluster order to cluster_row_order.csv and cluster_col_order.csv")
    finally:
        plt.close(g.fig)


def parse_args():
    p = argparse.ArgumentParser(description="Signature-Based AA Variant Modeling")
    p.add_argument('--signature_vector', required=True,
                        help='CSV: ID + 12 signature proportions')
    p.add_argument('--cosmic', action='store_true',
               help='Treat --signature_vector as a COSMIC 96-context table (Type + one signature column).')
    p.add_argument('--out_dir', required=True,
                        help='Directory for outputs')
    p.add_argument('--step', choices=['model','heatmap','both'], default='both',
                        help='Step to run')
    p.add_argument('--log_level', choices=['DEBUG','INFO','WARNING','ERROR'], default='INFO',
                        help='Logging level')
 
    # Signature comparison flags
    p.add_argument(
        "--compare_to_cosmic",
        help="Path to a COSMIC signature .txt file OR a folder of such files (e.g., COSMIC_v3.4_*_GRCh38.txt)."
    )
    p.add_argument(
        "--observed_csv",
        help="CSV with rows = samples, columns = mutation contexts (e.g., 96 SBS contexts like A[C>A]A)."
    )
    p.add_argument(
        "--heatmap_metric",
        default="cosine",
        choices=["cosine","pearson","l1","l2"],
        help="Similarity/distance metric for the heatmap (default: cosine)."
    )
    p.add_argument(
        "--heatmap_title",
        default="Observed vs COSMIC",
        help="Title for the heatmap."
    )
 
    # ---- plotting controls ----
    p.add_argument('--no-plot', action='store_true',
                        help='Skip heatmap even if step includes heatmap')
    p.add_argument('--cmap', default='viridis',
                        help='Matplotlib/Seaborn colormap name (e.g., viridis, magma, coolwarm)')
    p.add_argument('--figsize', default='12x8',
                        help='Figure size WxH (e.g., "12x8" or "12,8")')
    p.add_argument('--annot', action='store_true',
                        help='Annotate heatmap cells with values')
    p.add_argument('--linewidths', type=float, default=0.5,
                        help='Heatmap grid line width')
    p.add_argument('--linecolor', default='gray',
                        help='Heatmap grid line color')
    p.add_argument('--vmin', type=float, default=None,
                        help='Heatmap lower color bound')
    p.add_argument('--vmax', type=float, default=None,
                        help='Heatmap upper color bound')

    # clustering controls
    p.add_argument('--cluster', choices=['none','rows','cols','both'], default='both',
                   help='Hierarchical clustering for expected heatmap (default: both).')
    p.add_argument('--cluster_method', default='average',
                   choices=['average','complete','single','ward','weighted','centroid','median'],
                   help='Linkage method for clustering (default: average).')
    p.add_argument('--cluster_metric', default='euclidean',
                   choices=['euclidean','correlation','cosine','cityblock','chebyshev'],
                   help='Distance metric for clustering (default: euclidean).')
    p.add_argument('--cluster_standardize', default='none',
                   choices=['none','row','col'],
                   help='Standardize values across rows or cols before clustering (default: none).')


    # optional: change output CSV filename
    p.add_argument('--out_csv', default='expected_vectors.csv',
                        help='Filename for expected vectors CSV inside --out_dir')
    return p.parse_args()


def main():
    args = parse_args()
    init_logging(args.log_level)
    out = Path(args.out_dir)
    try:
        out.mkdir(parents=True, exist_ok=True)
    except PermissionError:
        logging.error(f"Permission denied: cannot create output directory {out}")
        return
        
    # Read signature vector
    try:
        if args.cosmic:
            # COSMIC: read raw (no index_col), then convert to 12-channel
            raw = pd.read_csv(args.signature_vector, sep=r"\t", engine="python", dtype=str)
            if "Type" in raw.columns and raw.shape[1] == 2:
                value_col = [c for c in raw.columns if c.lower() != "type"][0]
                s96 = pd.Series(pd.to_numeric(raw[value_col], errors="coerce").values,
                                index=raw["Type"].astype(str).str.strip(), dtype=float).dropna()
                twelve = cosmic96_to_12(s96)
                df = pd.DataFrame([twelve.values], index=["COSMIC"], columns=_TWELVE)
            else:
                raise ValueError("--signature_vector with --cosmic must be a two-column table (Type + one signature).")
        else:
            # 12-channel CSV: ID is the index
            df = pd.read_csv(args.signature_vector, index_col=0)
    except Exception as e:
        logging.error(f"Failed to read signature vector: {e}")
        return


    # Validate columns (must be the 12 contexts)
    try:
        validate_context_df(df)
    except ValueError as e:
        logging.error(e)
        return
    
    # If we need to build the model, do it and write expected vectors
    if args.step in ['model', 'both']:
        exp_df = build_expected(df)

        # Ensure there is a concrete 'ID' column exactly once
        if 'ID' in exp_df.columns:
            exp_df = exp_df.reset_index(drop=True)
        else:
            idx_name = exp_df.index.name  # could be 'ID', 'index', or None
            exp_df = exp_df.reset_index()
            exp_df = exp_df.rename(columns={(idx_name if idx_name else 'index'): 'ID'})

        csv_path = out / args.out_csv
        try:
            exp_df.to_csv(csv_path, index=False)
            logging.info(f"Wrote expected vectors -> {csv_path}")
        except PermissionError:
            logging.error(f"Permission denied: cannot write {csv_path}")
            return

        # If we need the heatmap, load (or reuse) expected vectors and plot
    if args.step in ['heatmap', 'both'] and not args.no_plot:
        try:
            if 'exp_df' not in locals():
                exp_df = pd.read_csv(out / args.out_csv)
            # Non-clustered (kept for convenience)
            plot_heatmap(
                exp_df.set_index('ID'), out,
                cmap=args.cmap,
                figsize=_parse_figsize(args.figsize),
                annot=args.annot,
                linewidths=args.linewidths,
                linecolor=args.linecolor,
                vmin=args.vmin,
                vmax=args.vmax,
            )
            # Clustered
            plot_heatmap_clustered(
                exp_df.set_index('ID'), out,
                cmap=args.cmap,
                figsize=_parse_figsize(args.figsize),
                vmin=args.vmin,
                vmax=args.vmax,
                cluster=args.cluster,
                method=args.cluster_method,
                metric=args.cluster_metric,
                standardize=args.cluster_standardize,
            )
        except FileNotFoundError:
            logging.error(f"{args.out_csv} not found; run with --step model or --step both first.")
        except Exception as e:
            logging.error(f"Heatmap failed: {e}")


    # --- Optional: compare observed vectors to COSMIC signatures and make a heatmap ---
    if args.compare_to_cosmic and args.observed_csv:
        try:
            # Load observed
            obs_path = Path(args.observed_csv)
            if not obs_path.exists():
                logging.error(f"[cosmic] observed CSV not found: {obs_path}")
            else:
                obs = pd.read_csv(obs_path, dtype=float)
                # If first column is IDs, use it as index
                if obs.columns[0].lower() in {"id","sample","sample_id","case","case_id"}:
                    obs.set_index(obs.columns[0], inplace=True)
                # Load signatures
                sigs = load_signatures_from_path(args.compare_to_cosmic)
                logging.info(f"[cosmic] Loaded {len(sigs)} signature(s) from {args.compare_to_cosmic}")
    
                out_png = Path(args.out_dir) / f"obs_vs_cosmic_{args.heatmap_metric}.png"
                out_csv = Path(args.out_dir) / f"obs_vs_cosmic_{args.heatmap_metric}.csv"
                compare_signatures_heatmap(
                    observed=obs,
                    signatures=sigs,
                    out_png=out_png,
                    out_csv=out_csv,
                    metric=args.heatmap_metric,
                    normalize_observed=True,
                    normalize_signatures=True,
                    title=args.heatmap_title,
                )
                logging.info(f"[cosmic] Wrote heatmap -> {out_png}")
                logging.info(f"[cosmic] Wrote matrix  -> {out_csv}")
        except Exception as e:
            logging.warning(f"[cosmic] comparison failed: {e}")
    elif args.compare_to_cosmic and not args.observed_csv:
        logging.warning("[cosmic] --compare_to_cosmic was provided but --observed_csv was not.")
         

if __name__=='__main__':
    main()
