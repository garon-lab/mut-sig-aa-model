#!/usr/bin/env python3
from __future__ import annotations
"""
MULTIOMIC INTEGRATION

This pipeline uses a manifest TSV to integrate multiple omics layers (DNA required; RNA, methylation, protein, and copy number optional)
for each case-id listed in a manifest. Can be used with an existing unifed manifest built with make_manifest.py or separate manifests. 
Output directory contains {case-id}_integrated.csv and out_dir/dna/{case-id}.csv. 
Note protein files must be preprocessed and in the format {case-id}.csv (see protein_preprocessor.py). 
DNA files similarly can be modified with the dna_preprocessor.py.

Dependencies: csv, argparse, shutil, pandas, numpy

Project layout (root = --folder {project}):
  {project}/
    dna/      {case_id}.csv (or similarly named)
    rna/      (subfolders per reference manifest) / files referenced by manifest
    ch3/      (subfolders per reference manifest) / files referenced by manifest
    cnv/      (subfolders per reference manifest) / files referenced by manifest
    protein/  {case_id}.csv (or similarly named)

Manifest (case list):
  - plain text: one case ID per line (comments '#' allowed), OR
  - CSV/TSV: first column is case_id

Steps:
  --step all|dna|rna|ch3|cnv|protein
  
  Running 'dna' will also chain RNA→CH3→Protein→CNV integration into DNA outputs,
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
    --folder <project root directory> \
    --case-id C3N-001 \
    --manifest <manifest file> \
    --out_dir <output directory> \
    --ref_zip <reference directory> \
    --step <all|dna|rna|ch3|cnv|protein>

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
    --emit_qc             Saves


Notes
Only the final protein files (with SNV/SNP, RNA, CH3, protein, and CNV) are written to: OUTPUT_DIR/{case-ID}.csv

"""
import argparse
import logging
import os
import re
import shutil
import tempfile
import gzip
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Optional, Tuple, List

import pandas as pd

logging.basicConfig(level=logging.INFO, format='[multiomic_integration] %(levelname)s %(message)s')

# ---------------------- helpers ----------------------

def _norm_ensg_ser(s: pd.Series) -> pd.Series:
    """Extract versionless Ensembl IDs (ENSG#########)."""
    return s.astype(str).str.extract(r'(ENSG\d+)', expand=False)

def _read_any_table(path: Path, sep: Optional[str] = None, names: Optional[List[str]] = None) -> pd.DataFrame:
    if sep is None:
        sep = "," if path.suffix.lower() == ".csv" else "\t"
    try:
        return pd.read_csv(path, sep=sep, engine="c", dtype=str, low_memory=False,
                           header=None if names else "infer", names=names)
    except Exception:
        return pd.read_csv(path, sep=sep, engine="python", dtype=str,
                           header=None if names else "infer", names=names)

def _safe_to_numeric(s: pd.Series) -> pd.Series:
    return pd.to_numeric(s, errors="coerce")

def _ensure_dir(p: Path):
    p.mkdir(parents=True, exist_ok=True)

def _mod_dir(folder: Path, primary: str, alt: Optional[str] = None) -> Path:
    """Prefer primary subdir if it exists, otherwise alt, otherwise primary."""
    p = folder / primary
    if p.exists():
        return p
    if alt:
        q = folder / alt
        if q.exists():
            return q
    return p

def _load_manifest_df(mod_dir: Path, manifest: Optional[str], default_name: str) -> Optional[pd.DataFrame]:
    """Load a case→file manifest. Accepts absolute path, or path relative to mod_dir.
    Normalizes a few common column names."""
    candidates: List[Path] = []
    if manifest:
        mp = Path(manifest)
        candidates += [mp, mod_dir / mp.name, mod_dir / manifest]
    candidates.append(mod_dir / default_name)

    man_path = next((p for p in candidates if p and p.exists()), None)
    if man_path is None:
        return None

    df = _read_any_table(man_path, sep="\t")
    # normalize column names
    norm = {re.sub(r"\s+", " ", c.lower().replace("-", " ").replace("_", " ").strip()): c for c in df.columns}
    def pick(*keys):
        for k in keys:
            if k in norm: return norm[k]
        return None

    case_c = pick("case id", "case", "sample id", "sample")
    fileid_c = pick("file id", "id")
    fname_c  = pick("file name", "filename", "name")
    path_c   = pick("path", "file path", "filepath", "relpath", "full path")

    if not case_c or not (fileid_c or fname_c or path_c):
        logging.warning(f"Manifest {man_path} missing required columns (need Case ID and one of File ID / File Name / Path).")
        return None

    out = pd.DataFrame({"Case ID": df[case_c].astype(str).str.strip()})
    if fileid_c: out["File ID"] = df[fileid_c].astype(str).str.strip()
    if fname_c:  out["File Name"] = df[fname_c].astype(str).str.strip()
    if path_c:   out["Path"] = df[path_c].astype(str).str.strip()
    return out.dropna(subset=["Case ID"])

def _resolve_data_root(mod_dir: Path, row: pd.Series) -> Optional[Path]:
    """Return a directory to search OR a direct file path if Path points to the file."""
    # Explicit path wins
    p = Path(str(row.get("Path"))) if "Path" in row and pd.notna(row.get("Path")) else None
    if p:
        if not p.is_absolute():
            p = (mod_dir / p)
        if p.exists():
            return p

    # GDC-style by File ID (a folder)
    if "File ID" in row and pd.notna(row.get("File ID")):
        candidate = mod_dir / str(row["File ID"])
        if candidate.exists():
            return candidate

    # Fallback: search by file name (returns the file path)
    if "File Name" in row and pd.notna(row.get("File Name")):
        for hit in mod_dir.rglob(str(row["File Name"])):
            return hit  # file path
    return None


# ---------------------- reference (maps only) -----------------------

EXPECTED_BASENAMES = {
    "esng_gene-sym.txt",
    "gene-sym_name.txt",
    "uni-ensg_all.txt",
    "uni-np_all.txt",
}

def _pick_best_subdir(tmp_root: Path) -> Path:
    subs = [d for d in tmp_root.iterdir() if d.is_dir()]
    for d in subs:
        if d.name.lower() == "reference":
            return d
    for d in subs:
        names = {p.name for p in d.iterdir() if p.is_file()}
        if EXPECTED_BASENAMES & names:
            return d
    return subs[0] if subs else tmp_root

def resolve_reference_dir(ref_dir: Optional[str], ref_zip: Optional[str]) -> Optional[Path]:
    if ref_dir:
        p = Path(ref_dir)
        if not p.exists():
            raise FileNotFoundError(f"--ref_dir does not exist: {p}")
        return p
    if ref_zip:
        zp = Path(ref_zip)
        if not zp.exists():
            raise FileNotFoundError(f"--ref_zip not found: {zp}")
        if zp.suffix.lower() != ".zip":
            raise ValueError(f"--ref_zip must be a .zip file: {zp}")
        tmp = Path(tempfile.mkdtemp(prefix="reference_"))
        shutil.unpack_archive(str(zp), str(tmp), format="zip")
        best = _pick_best_subdir(tmp)
        logging.info(f"Extracted reference to: {best}")
        return best
    return None

# ---------------------------- manifest ------------------------------

def load_case_ids(manifest_path: str) -> List[str]:
    s = pd.read_table(manifest_path, usecols=[0], header=0, dtype=str,
                      na_values=["", "NA", "NaN"]).iloc[:, 0]
    s = s.dropna().astype(str).str.strip()
    header_like = {"case-id","case_id","id","sample-id","sample_id","Case-ID","Case_ID","ID","Sample-ID","Sample_ID"}
    s = s[~s.str.lower().isin({x.lower() for x in header_like})]
    return [x for x in s.unique() if x]

def _read_gdc_manifest(mod_dir: Path, basename: str) -> Optional[pd.DataFrame]:
    m = mod_dir / basename
    if not m.exists():
        return None
    df = _read_any_table(m, sep="\t")
    lc = {c.lower().replace("-", " ").replace("_", " "): c for c in df.columns}
    case_col = lc.get("case id")
    file_col = lc.get("file id")
    name_col = lc.get("file name")
    if case_col and file_col:
        keep = [case_col, file_col] + ([name_col] if name_col else [])
        return df[keep].rename(columns={case_col: "Case ID", file_col: "File ID"})
    return None


# ----------------------------- DNA ---------------------------------

def _add_ENSG_columns(df: pd.DataFrame, case_id: str) -> pd.DataFrame:
    """Ensure both ENSGene (full) and ENSGene_core exist; insert ENSGene right after INFO if INFO exists."""
    df = df.copy()
    if "INFO" in df.columns:
        info_clean = df["INFO"].astype(str).str.strip().str.strip('"').str.strip("'")
        ensg_full = info_clean.str.extract(r'(ENSG\d+(?:\.\d+)?)', expand=False)
        if "ENSGene" in df.columns:
            df["ENSGene"] = df["ENSGene"].where(df["ENSGene"].notna(), ensg_full)
        else:
            info_idx = df.columns.get_loc("INFO")
            df.insert(info_idx + 1, "ENSGene", ensg_full)
        logging.info(f"[DNA] {case_id}: ENSG extracted for {int(df['ENSGene'].notna().sum())} / {len(df)} rows")
    else:
        if "ENSGene" not in df.columns:
            df.insert(len(df.columns), "ENSGene", pd.Series(index=df.index, dtype=object))
            logging.warning(f"[DNA] {case_id}: INFO missing; added empty ENSGene column.")
    # Always build core
    if "ENSGene" in df.columns:
        df["ENSGene_core"] = _norm_ensg_ser(df["ENSGene"])
    elif "INFO" in df.columns:
        df["ENSGene_core"] = _norm_ensg_ser(df["INFO"])
    else:
        df["ENSGene_core"] = pd.Series(index=df.index, dtype=object)
    return df

def _try_per_case_file(folder: Path, modality: str, case_id: str) -> Optional[Path]:
    mod_dir = folder / modality
    for ext in (".csv", ".tsv", ".txt"):
        p = mod_dir / f"{case_id}{ext}"
        if p.exists():
            return p
    return None

def load_and_enhance_dna(case_id: str, folder: Path, out_dir: Path, ref_dir: Optional[Path], dna_rows: str) -> Path:
    inp = _try_per_case_file(folder, "dna", case_id)
    if inp is None:
        tried = [str(folder / "dna" / f"{case_id}{ext}") for ext in (".csv",".tsv",".txt")]
        raise FileNotFoundError(f"DNA file not found for {case_id}. Tried: {', '.join(tried)}")
    df = _read_any_table(inp)
    df = _add_ENSG_columns(df, case_id)
    if dna_rows == "ensg_only":
        before = len(df)
        df2 = df[df["ENSGene_core"].notna()].copy()
        if df2.empty:
            logging.warning(f"[DNA] {case_id}: ENSG-only filter removed all rows; keeping all rows.")
        else:
            df = df2
        logging.info(f"[DNA] {case_id}: filtered to ENSG-only rows: {len(df)} / {before}")
    df = integrate_gene_annotations(df, ref_dir, case_id)
    df = integrate_uniprot_np(df, ref_dir, case_id)
    out_dna = out_dir / "dna"; _ensure_dir(out_dna)
    out_fp = out_dna / f"{case_id}.csv"
    df.to_csv(out_fp, index=False)
    return out_fp

# ----------------------- modality I/O (GDC-style) -------------------

def _find_in_subdir_by_patterns(root: Path, patterns: List[str]) -> Optional[Path]:
    if not root.exists(): return None
    for p in root.rglob("*"):
        if p.is_file():
            name = p.name.lower()
            if all(pat in name for pat in patterns):
                return p
    return None

def _read_first_two_col_table(path: Path, gz_ok: bool = True) -> pd.DataFrame:
    if path.suffix == ".gz" and gz_ok:
        with gzip.open(path, "rt") as fh:
            return pd.read_csv(fh, sep="\t", header=None, names=["col1","col2"], dtype=str)
    else:
        return pd.read_csv(path, sep="\t", header=None, names=["col1","col2"], dtype=str)

# -------------------------- reference maps --------------------------

def _find_ref_file(ref_dir: Optional[Path], basename: str) -> Optional[Path]:
    if ref_dir is None: return None
    direct = ref_dir / basename
    if direct.exists(): return direct
    for p in ref_dir.rglob(basename):
        if "__MACOSX" in p.parts: continue
        return p
    return None

def integrate_gene_annotations(base_df: pd.DataFrame, ref_dir: Optional[Path], case_id: str) -> pd.DataFrame:
    if ref_dir is None: return base_df
    df = base_df.copy()
    sym_fp = _find_ref_file(ref_dir, "esng_gene-sym.txt")
    name_fp = _find_ref_file(ref_dir, "gene-sym_name.txt")
    if sym_fp:
        m = _read_any_table(sym_fp, sep="\t", names=["Gene","Ensembl"]).drop_duplicates()
        m["Ensembl_core"] = _norm_ensg_ser(m["Ensembl"])
        m = m.drop_duplicates(subset=["Ensembl_core"], keep="first")
        df = df.merge(m[["Gene","Ensembl_core"]], how="left",
                      left_on="ENSGene_core", right_on="Ensembl_core")
    else:
        logging.warning(f"[DNA] {case_id}: Gene symbol map not found in reference")
    if name_fp and "Gene" in df.columns:
        n = _read_any_table(name_fp, sep="\t", names=["symbol","name"]).drop_duplicates()
        n = n.drop_duplicates(subset=["symbol"], keep="first")
        df = df.merge(n, how="left", left_on="Gene", right_on="symbol")
    return df

def integrate_uniprot_np(base_df: pd.DataFrame, ref_dir: Optional[Path], case_id: str) -> pd.DataFrame:
    if ref_dir is None: return base_df
    df = base_df.copy()
    ue_fp = _find_ref_file(ref_dir, "uni-ensg_all.txt")
    unp_fp = _find_ref_file(ref_dir, "uni-np_all.txt")
    if ue_fp:
        ue = _read_any_table(ue_fp, sep="\t", names=["From","To"]).drop_duplicates()
        ue["To_core"] = _norm_ensg_ser(ue["To"])
        ue = ue.drop_duplicates(subset=["To_core"], keep="first")
        df = df.merge(ue[["From","To_core"]], how="left",
                      left_on="ENSGene_core", right_on="To_core")
    else:
        logging.warning(f"[DNA] {case_id}: UniProt→ENSG map not found")
    if unp_fp and "From" in df.columns:
        unp = _read_any_table(unp_fp, sep="\t", names=["From","NP"]).drop_duplicates()
        unp = unp.drop_duplicates(subset=["From"], keep="first")
        df["From"] = df["From"].astype(str).str.strip()
        unp["From"] = unp["From"].astype(str).str.strip()
        df = df.merge(unp, how="left", on="From")
    elif not unp_fp:
        logging.warning(f"[DNA] {case_id}: UniProt→NP map not found")
    return df

# ----------------------- modality merges ----------------------------

def _rna_from_gdc(folder: Path, case_id: str, rna_manifest: Optional[str]) -> Optional[pd.DataFrame]:
    mod_dir = folder / "rna"
    man = _load_manifest_df(mod_dir, rna_manifest, "gdc-rna.tsv")
    if man is None: return None
    rows = man[man["Case ID"] == case_id]
    if rows.empty: return None

    root_or_file = _resolve_data_root(mod_dir, rows.iloc[0])
    if root_or_file is None:
        logging.warning(f"[RNA] {case_id}: cannot resolve path from manifest.")
        return None

    f = root_or_file if root_or_file.is_file() else _find_in_subdir_by_patterns(root_or_file, ["htseq", "counts"])
    if f is None or not f.exists():
        logging.warning(f"[RNA] {case_id}: no htseq counts file under {root_or_file}")
        return None

    df = _read_first_two_col_table(f, gz_ok=True).rename(columns={"col1":"Ensembl_full","col2":"Count"})
    df["Ensembl_core"] = _norm_ensg_ser(df["Ensembl_full"])
    df = df.dropna(subset=["Ensembl_full"]).drop_duplicates(subset=["Ensembl_full"], keep="first")
    df["Count"] = _safe_to_numeric(df["Count"])
    return df


def _cn_from_gdc(folder: Path, case_id: str, cn_manifest: Optional[str]) -> Optional[pd.DataFrame]:
    mod_dir = _mod_dir(folder, "cnv", alt="cn")
    man = _load_manifest_df(mod_dir, cn_manifest, "gdc-cn.tsv")
    if man is None: return None
    rows = man[man["Case ID"] == case_id]
    if rows.empty: return None

    root_or_file = _resolve_data_root(mod_dir, rows.iloc[0])
    if root_or_file is None:
        logging.warning(f"[CNV] {case_id}: cannot resolve path from manifest.")
        return None

    f = root_or_file if root_or_file.is_file() else _find_in_subdir_by_patterns(root_or_file, ["gene_level", "copy_number_variation"])
    if f is None or not f.exists():
        logging.warning(f"[CNV] {case_id}: no gene_level.copy_number_variation under {root_or_file}")
        return None

    df = _read_any_table(f, sep="\t")
    lc = {c.lower().replace("-", " ").replace("_", " "): c for c in df.columns}
    gid = lc.get("gene id"); cn = lc.get("copy number")
    if gid is None or cn is None:
        logging.warning(f"[CNV] {case_id}: expected gene_id and copy_number in {f.name}")
        return None
    out = df[[gid, cn]].rename(columns={gid:"Ensembl_full", cn:"copy_number"})
    out["Ensembl_core"] = _norm_ensg_ser(out["Ensembl_full"])
    out = out.dropna(subset=["Ensembl_full"]).drop_duplicates(subset=["Ensembl_full"], keep="first")
    out["copy_number"] = _safe_to_numeric(out["copy_number"])
    return out

def _ch3_from_gdc(folder: Path, case_id: str, ch3_manifest: Optional[str]) -> Optional[pd.DataFrame]:
    mod_dir = folder / "ch3"
    man = _load_manifest_df(mod_dir, ch3_manifest, "gdc-ch3.tsv")
    if man is None: return None
    rows = man[man["Case ID"] == case_id]
    if rows.empty: return None

    root_or_file = _resolve_data_root(mod_dir, rows.iloc[0])
    if root_or_file is None:
        logging.warning(f"[CH3] {case_id}: cannot resolve path from manifest.")
        return None

    # If the manifest points directly to a table, use it; else search for a cg/beta TSV
    candidate = None
    if root_or_file.is_file():
        candidate = root_or_file
    else:
        for p in root_or_file.rglob("*.tsv"):
            try:
                head = pd.read_csv(p, sep="\t", dtype=str, nrows=16)
                lc = {c.lower(): c for c in head.columns}
                if (lc.get("cg") or lc.get("ilmnid") or lc.get("composite element ref")) and \
                   (lc.get("beta") or lc.get("beta_value") or lc.get("beta value")):
                    candidate = p; break
            except Exception:
                continue
    if candidate is None:
        logging.warning(f"[CH3] {case_id}: could not locate cg/beta table under {root_or_file}")
        return None

    full = pd.read_csv(candidate, sep="\t", dtype=str)
    lc = {c.lower(): c for c in full.columns}
    probe_col = lc.get("cg") or lc.get("ilmnid") or lc.get("composite element ref")
    beta_col  = lc.get("beta") or lc.get("beta_value") or lc.get("beta value")
    out = full[[probe_col, beta_col]].rename(columns={probe_col:"probe", beta_col:"beta"})
    out["beta"] = _safe_to_numeric(out["beta"])
    return out.dropna(subset=["probe"])


def _load_ch3_map(ch3_map: Optional[str], ref_dir: Optional[Path],
                  ch3_probe_col: Optional[str], ch3_ensg_col: Optional[str], ch3_symbol_col: Optional[str]) -> Optional[pd.DataFrame]:
    path: Optional[Path] = None
    if ch3_map:
        path = Path(ch3_map)
    else:
        for name in ["probe-ensg.tsv","probe_ensg.tsv","illumina_probes_ensg.tsv","ch3_probe_map.tsv","ch3_map.tsv"]:
            guess = _find_ref_file(ref_dir, name) if ref_dir else None
            if guess and guess.exists(): path = guess; break
    if path is None or not path.exists():
        return None
    df = _read_any_table(path)
    lc = {c.lower(): c for c in df.columns}
    probe_c = ch3_probe_col or lc.get("probe") or lc.get("cg") or lc.get("ilmnid") or lc.get("composite element ref")
    ensg_c  = ch3_ensg_col  or lc.get("ensg") or lc.get("ensembl") or lc.get("ensembl_gene_id") or lc.get("gene_id")
    sym_c   = ch3_symbol_col or lc.get("gene") or lc.get("symbol") or lc.get("gene_symbol")
    if not probe_c or not (ensg_c or sym_c):
        logging.warning(f"[CH3] Map {path} missing probe+gene columns; specify --ch3_probe_col/--ch3_ensg_col/--ch3_symbol_col.")
        return None
    m = df.rename(columns={probe_c:"probe"}).copy()
    if ensg_c:
        m["ENSGene_core"] = _norm_ensg_ser(m[ensg_c])
    else:
        sym_fp = _find_ref_file(ref_dir, "esng_gene-sym.txt") if ref_dir else None
        if sym_fp is None:
            logging.warning("[CH3] No symbol map (esng_gene-sym.txt) in reference; cannot derive ENSG from symbol.")
            return None
        sym = _read_any_table(sym_fp, sep="\t", names=["Gene","Ensembl"]).drop_duplicates()
        sym["Ensembl_core"] = _norm_ensg_ser(sym["Ensembl"])
        sym = sym.drop_duplicates(subset=["Gene"], keep="first")
        m = m.rename(columns={sym_c:"Gene"}).merge(
            sym[["Gene","Ensembl_core"]], how="left", on="Gene"
        ).rename(columns={"Ensembl_core":"ENSGene_core"})
    m = m.dropna(subset=["probe","ENSGene_core"]).drop_duplicates(subset=["probe"], keep="first")
    return m[["probe","ENSGene_core"]]

def integrate_rna(base_df: pd.DataFrame, folder: Path, case_id: str,
                  ensg_join_mode: str, synth_overlap: bool, rna_manifest: Optional[str]) -> pd.DataFrame:
    df = base_df.copy()
    rna = _rna_from_gdc(folder, case_id, rna_manifest)
    if rna is None or rna.empty: return df
    left_key  = "ENSGene_core" if ensg_join_mode == "core" else "ENSGene"
    right_key = "Ensembl_core" if ensg_join_mode == "core" else "Ensembl_full"
    if synth_overlap:
        need = set(df[left_key].dropna().unique()) - set(rna[right_key].dropna().unique())
        if need:
            add = pd.DataFrame({
                "Ensembl_full": [x if ensg_join_mode=="exact" else (f"{x}.1" if not re.search(r'\.\d+$', str(x)) else str(x)) for x in need],
                "Ensembl_core": list(need), "Count": 100
            })
            rna = pd.concat([rna, add], ignore_index=True)
    df = df.merge(rna[[right_key,"Count"]], how="left", left_on=left_key, right_on=right_key)
    logging.info(f"[RNA] {case_id}: matched Count for {int(df['Count'].notna().sum())} rows (join={ensg_join_mode})")
    return df

def integrate_cnv(base_df: pd.DataFrame, folder: Path, case_id: str,
                  ensg_join_mode: str, synth_overlap: bool, cn_manifest: Optional[str]) -> pd.DataFrame:
    df = base_df.copy()
    cn = _cn_from_gdc(folder, case_id, cn_manifest)
    if cn is None or cn.empty: return df
    left_key  = "ENSGene_core" if ensg_join_mode == "core" else "ENSGene"
    right_key = "Ensembl_core" if ensg_join_mode == "core" else "Ensembl_full"
    if synth_overlap:
        need = set(df[left_key].dropna().unique()) - set(cn[right_key].dropna().unique())
        if need:
            add = pd.DataFrame({
                "Ensembl_full": [x if ensg_join_mode=="exact" else (f"{x}.1" if not re.search(r'\.\d+$', str(x)) else str(x)) for x in need],
                "Ensembl_core": list(need), "copy_number": 2
            })
            cn = pd.concat([cn, add], ignore_index=True)
    df = df.merge(cn[[right_key,"copy_number"]], how="left", left_on=left_key, right_on=right_key)
    logging.info(f"[CNV] {case_id}: matched copy_number for {int(df['copy_number'].notna().sum())} rows (join={ensg_join_mode})")
    return df

def integrate_ch3(base_df: pd.DataFrame, folder: Path, case_id: str, ensg_join_mode: str,
                  ch3_map: Optional[str], ref_dir: Optional[Path], ch3_probe_col: Optional[str],
                  ch3_ensg_col: Optional[str], ch3_symbol_col: Optional[str], ch3_agg: str, ch3_manifest: Optional[str]) -> pd.DataFrame:
    df = base_df.copy()
    ch3 = _ch3_from_gdc(folder, case_id, ch3_manifest)
    if ch3 is None or ch3.empty: return df
    m = _load_ch3_map(ch3_map, ref_dir, ch3_probe_col, ch3_ensg_col, ch3_symbol_col)
    if m is None or m.empty:
        logging.warning(f"[CH3] {case_id}: No usable probe→ENSG map; skipping CH3.")
        return df
    ch3m = ch3.merge(m, how="left", on="probe").dropna(subset=["ENSGene_core"])
    if ch3m.empty:
        logging.warning(f"[CH3] {case_id}: probe map didn’t cover this case’s probes; skipping.")
        return df
    agg = (ch3m.groupby("ENSGene_core", as_index=False)["beta"].median() if ch3_agg=="median"
           else ch3m.groupby("ENSGene_core", as_index=False)["beta"].mean())
    agg = agg.rename(columns={"beta":"beta_val"})
    left_key = "ENSGene_core" if ensg_join_mode == "core" else "ENSGene"
    df = df.merge(agg, how="left", left_on=left_key, right_on="ENSGene_core")
    logging.info(f"[CH3] {case_id}: matched beta_val for {int(df['beta_val'].notna().sum())} rows (join={ensg_join_mode})")
    return df

def integrate_protein(base_df: pd.DataFrame, folder: Path, case_id: str, synth_overlap: bool) -> pd.DataFrame:
    p = _try_per_case_file(folder, "protein", case_id)
    if p is None: return base_df
    prot = _read_any_table(p, sep=",").drop_duplicates()
    if "NP" not in prot.columns:
        logging.warning(f"[PROTEIN] {case_id}: protein file lacks NP; skipping.")
        return base_df
    df = base_df.copy()
    base_np = set(df.get("NP", pd.Series([], dtype=str)).dropna().astype(str).str.strip().unique().tolist())
    prot_np = set(prot["NP"].dropna().astype(str).str.strip().unique().tolist())
    need = base_np - prot_np
    if synth_overlap and need:
        logging.info(f"[PROTEIN] {case_id}: synthesizing {len(need)} NP rows")
        prot = pd.concat([prot, pd.DataFrame({"NP": list(need), "SEQ":"TESTPEPTIDE","EV":"1e-5","INT":"1000/0.0"})], ignore_index=True)
    prot = prot.drop_duplicates(subset=["NP"], keep="first")
    df["NP"] = df.get("NP", pd.Series(index=df.index, dtype=object)).astype(str).str.strip()
    prot["NP"] = prot["NP"].astype(str).str.strip()
    df = df.merge(prot, how="left", on="NP", suffixes=("", "_protein"))
    logging.info(f"[PROTEIN] {case_id}: matched NP for {int(df['SEQ'].notna().sum()) if 'SEQ' in df.columns else 0} rows")
    return df

# -------------------------- dedup & QC -------------------------------

def _final_dedup(df: pd.DataFrame, level: str, key_cols: Optional[List[str]]) -> pd.DataFrame:
    if level == "none": return df
    if level == "row":  return df.drop_duplicates()
    if level == "key" and key_cols:
        cols = [c for c in key_cols if c in df.columns]
        if cols: return df.drop_duplicates(subset=cols, keep="first")
    return df

# -------------------------- per-case worker -------------------------

def process_one_case(case_id: str, folder: str, out_dir: str, ref_dir: Optional[str],
                     step: str, dna_rows: str, dedup_level: str, dedup_key_cols: Optional[List[str]],
                     emit_qc: bool, ensg_join_mode: str, synth_overlap: bool,
                     ch3_map: Optional[str], ch3_probe_col: Optional[str], ch3_ensg_col: Optional[str],
                     ch3_symbol_col: Optional[str], ch3_agg: str,
                     skip_rna: bool, skip_ch3: bool, skip_protein: bool, skip_cn: bool,
                     rna_manifest: Optional[str], ch3_manifest: Optional[str], cn_manifest: Optional[str]
                     ) -> Tuple[str, Optional[str]]:
    try:
        folder_p = Path(folder); out_dir_p = Path(out_dir); ref_p = Path(ref_dir) if ref_dir else None

        # Prepare DNA if requested
        if step in ("dna","all"):
            load_and_enhance_dna(case_id, folder_p, out_dir_p, ref_p, dna_rows)

        # Always start merges from the DNA output table
        base = _read_any_table(out_dir_p / "dna" / f"{case_id}.csv")

        # Decide which integrations to run
        chain = (step == "dna")
        do_rna = ((step in ("rna","all")) or chain) and (not skip_rna)
        do_ch3 = ((step in ("ch3","all")) or chain) and (not skip_ch3)
        do_pro = ((step in ("protein","all")) or chain) and (not skip_protein)
        do_cnv = ((step in ("cnv","all")) or chain) and (not skip_cn)

        if do_rna:
            base = integrate_rna(base, folder_p, case_id, ensg_join_mode, synth_overlap, rna_manifest)
          
        if do_ch3:
            base = integrate_ch3(base, folder_p, case_id, ensg_join_mode, ch3_map, ref_p,
                                 ch3_probe_col, ch3_ensg_col, ch3_symbol_col, ch3_agg, ch3_manifest)
        if do_pro:
            base = integrate_protein(base, folder_p, case_id, synth_overlap)
          
        if do_cnv:
            base = integrate_cnv(base, folder_p, case_id, ensg_join_mode, synth_overlap, cn_manifest)


        base = _final_dedup(base, dedup_level, dedup_key_cols)

        out_fp = Path(out_dir) / f"{case_id}_integrated.csv"
        base.to_csv(out_fp, index=False)

        if emit_qc:
            qc = {
                "case_id": case_id,
                "n_rows_out": len(base),
                "rna_nonnull": int(base["Count"].notna().sum()) if "Count" in base.columns else 0,
                "cnv_nonnull": int(base["copy_number"].notna().sum()) if "copy_number" in base.columns else 0,
                "ch3_nonnull": int(base["beta_val"].notna().sum()) if "beta_val" in base.columns else 0,
                "prot_nonnull": int(base["SEQ"].notna().sum()) if "SEQ" in base.columns else 0,
            }
            qcd = Path(out_dir) / "qc"; _ensure_dir(qcd)
            pd.DataFrame([qc]).to_csv(qcd / f"{case_id}_qc.tsv", sep="\t", index=False)

        return (case_id, None)
    except Exception as e:
        return (case_id, str(e))


# ------------------------------- main -------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Multiomic integration (11k3 fixed): core ENSG join, maps, CH3, synthesis.")
    p.add_argument("--folder", required=True)
    p.add_argument("--manifest", required=True)
    p.add_argument("--out_dir", required=True)
    p.add_argument("--ref_dir", default=None)
    p.add_argument("--ref_zip", default=None)
    p.add_argument("--step", default="all", choices=["all","dna","rna","ch3","cnv","protein"])
    p.add_argument("--jobs", type=int, default=os.cpu_count() or 2)
    p.add_argument("--case-id", default=None, help="Process a single case-id instead of a manifest")

    # DNA handling / output controls
    p.add_argument("--dna_rows", default="ensg_only", choices=["ensg_only","all"])
    p.add_argument("--dedup_level", default="row", choices=["none","row","key"])
    p.add_argument("--dedup_key", default="CHROM,POS,REF,ALT,ENSGene,ENSGene_core")
    p.add_argument("--emit_qc", action="store_true")

    # Join behavior
    p.add_argument("--ensg_join_mode", default="core", choices=["core","exact"])
    p.add_argument("--synthesize_overlap_for_tests", action="store_true")

    # CH3 options
    p.add_argument("--ch3_map", default=None)
    p.add_argument("--ch3_probe_col", default=None)
    p.add_argument("--ch3_ensg_col", default=None)
    p.add_argument("--ch3_symbol_col", default=None)
    p.add_argument("--ch3_agg", default="mean", choices=["mean","median"])

    # Skip flags
    p.add_argument("--skip-rna", action="store_true", help="Skip RNA integration when chaining from --step dna or all.")
    p.add_argument("--skip-ch3", action="store_true", help="Skip CH3 integration when chaining from --step dna or all.")
    p.add_argument("--skip-protein", action="store_true", help="Skip protein integration when chaining from --step dna or all.")
    p.add_argument("--skip-cn", action="store_true", help="Skip CNV integration when chaining from --step dna or all.")

    # Manifests for raw data
    p.add_argument("--rna-manifest", dest="rna_manifest", default=None,
                   help="Path to RNA manifest (case→file). Defaults to <folder>/rna/gdc-rna.tsv")
    p.add_argument("--ch3-manifest", dest="ch3_manifest", default=None,
                   help="Path to CH3 manifest (case→file). Defaults to <folder>/ch3/gdc-ch3.tsv")
    p.add_argument("--cn-manifest",  dest="cn_manifest",  default=None,
                   help="Path to CNV manifest (case→file). Defaults to <folder>/(cnv|cn)/gdc-cn.tsv")


    return p.parse_args()


def main() -> None:
    args = parse_args()
    out_dir = Path(args.out_dir); _ensure_dir(out_dir)
    ref_dir_p = resolve_reference_dir(args.ref_dir, args.ref_zip)
    if ref_dir_p is not None:
        logging.info(f"Using reference dir for maps: {ref_dir_p}")
    case_ids = load_case_ids(args.manifest)
    dedup_key_cols = [c.strip() for c in args.dedup_key.split(",") if c.strip()] if args.dedup_level == "key" else None

    if args.case_id:
        case_ids = [args.case_id]
    else:
        case_ids = load_case_ids(args.manifest)

    logging.info(f"N = {len(case_ids)} samples")
    logging.info(f"Running step(s): {args.step}")
    with ProcessPoolExecutor(max_workers=args.jobs) as ex:
        futs = [ex.submit(
            process_one_case, cid, args.folder, args.out_dir,
            str(ref_dir_p) if ref_dir_p else None,
            args.step, args.dna_rows, args.dedup_level, dedup_key_cols, args.emit_qc,
            args.ensg_join_mode, args.synthesize_overlap_for_tests,
            args.ch3_map, args.ch3_probe_col, args.ch3_ensg_col, args.ch3_symbol_col, args.ch3_agg,
            args.skip_rna, args.skip_ch3, args.skip_protein, args.skip_cn,
            args.rna_manifest, args.ch3_manifest, args.cn_manifest
        ) for cid in case_ids]

        for cid in case_ids:
            logging.info(f"Processing {cid} [step={args.step}]")
        for f in as_completed(futs):
            cid, err = f.result()
            if err:
                logging.error(f"[{cid}]: {err}")
    logging.info("Done.")

if __name__ == "__main__":
    main()
