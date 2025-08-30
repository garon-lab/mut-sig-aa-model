#!/usr/bin/env python3
from __future__ import annotations
"""
MULTIOMIC INTEGRATION

This pipeline uses a manifest TSV to integrate multiple omics layers (DNA required; RNA, methylation, protein, and copy number optional)
for each case-id listed in a manifest. It produces a DNA-annotated table and an integrated table per case, with summary metadata for the run.

Features: 
1. DNA-first workflow loads per-case DNA CSVs, annotates with ENSG extraction, then optionally joins RNA, CH3, CNV, and protein.
2. Parallel execution via ProcessPoolExecutor (--jobs).
3. Flexible manifests that use a unified case list plus per-modality manifests (auto-detected by default) or point to explicit manifest paths.
4. Normal tissue support (optional) to incorporate multiple omics layers from matched normal tissue (--add-normal-dir) . 
5. Robust omics integration for CPTAC tables.
6. Clean metadata outputs - a run log and QC summary are created in `metadata/`.

Output directory contains {case-id}_integrated.csv and out_dir/dna/{case-id}.csv. 
Note protein files must be preprocessed and in the format {case-id}.csv (see protein_preprocessor.py). 
DNA files similarly can be modified with the dna_preprocessor.py.

Dependencies: argparse, csv, gzip, logging, os, re, shutil, sys, tempfile, concurrent.futures, pathlib, pandas, numpy

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

Per-modality manifests (optional; auto-detected if omitted):
  - `--rna-manifest`  (case → RNA file path or File ID)
  - `--ch3-manifest`  (case → CH3 file path or File ID)
  - `--cn-manifest`   (case → CNV file path or File ID)

Reference:
  - By default, we look for 'reference.zip' alongside the script or in CWD.
  - You can pass --ref_zip /path/to/reference.zip OR --ref_dir /path/to/extracted/reference.
  - Inside the reference, we expect manifest files (filenames configurable via flags):
      * --rna-manifest (default: rna_manifest.tsv)
      * --ch3-manifest (default: ch3_manifest.tsv)
      * --cn-manifest  (default: cn_manifest.tsv)

Steps:
  --step all|dna|rna|ch3|cnv|protein
  
  Running 'dna' will also chain RNA→CH3→Protein→CNV integration into DNA outputs,
  unless you pass skip flags (e.g., --skip-rna, --skip-ch3, --skip-protein, --skip-cn)

Join behavior:
- `--ensg_join_mode {core|exact}` (default: `core`):
  - `core`: join by versionless `ENSG#########` (recommended),
  - `exact`: require versioned exact matches (e.g., `ENSG00000123456.13`).
- CNV fast-path: looks for `gene_id` (extract ENSG core) and `copy_number`
  directly. If only gene symbols are present, a symbol→ENSG map from the
  reference is used.
- Advanced: `--synthesize_overlap_for_tests` can fill missing keys for testing
  (e.g., default CN=2).

Usage: 
Recommended
python multiomic_integration.py \
    --in_dir <project root directory> \
    --manifest <manifest file> \
    --out_dir <output directory> \
    --ref_zip <reference.zip> \
    
Optional Single File
python multiomic_integration.py \
    --in_dir <project root directory> \
    --case-id C3N-001 \
    --manifest <manifest file> \
    --out_dir <output directory> \
    --ref_dir <reference directory> \
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
        --step <all|dna|rna|ch3|cnv|protein>

Arguments:
Required
    --in_dir              Project root directory
    --manifest            Case list file - formats supported: plain text (one case-ID per line), CSV/TSV (case-IDs in first column)
    --out_dir             Directory to write final integrated files

Reference Options
    --ref_dir             Directory with references for multiomic mapping if not using included reference.zip
    --ref_zip             Path to a reference.zip archive. By default, script looks for one in current directory --ref_dir <path>

Control & Performance
  --step                   Pipeline step(s) to run (default: all)
  --jobs                   Parallel workers (default: os.cpu_count())
  --case-id                Process a single case instead of all in --manifest
  --ensg_join_mode         {core|exact} (default: core)
  --keep-join-cols         Keep intermediate join keys in outputs (default: drop)

Input Overrides
  --input_dna_dir          (default: <in_dir>/dna)
  --input_rna_dir          (default: <in_dir>/rna)
  --input_ch3_dir          (default: <in_dir>/ch3)
  --input_cn_dir           (default: <in_dir>/cnv or <in_dir>/cn)
  --input_protein_dir      (default: <in_dir>/protein)
  --rna-manifest           case→RNA lookup (auto-detected if omitted)
  --ch3-manifest           case→CH3 lookup (auto-detected if omitted)
  --cn-manifest            case→CNV lookup (auto-detected if omitted)

Normal
  --add-normal-dir, --normal-dir   Root folder containing matched normal layers

QC & Logging
  --emit_qc                (This build writes **summary + cohort QC** in metadata/
                           and does **not** emit per-case QC files.)
  --log-file               Custom log path (default: <out_dir>/metadata/multiomic_run.log)
  --summary-out            Custom summary path (default: <out_dir>/metadata/qc_summary.tsv)
  --summary-format         {tsv|csv} (default inferred from --summary-out)
  --qc-cohort-out          Custom cohort QC CSV (default: <out_dir>/metadata/qc_cohort.csv)

    --step                Step(s) to run [all|dna|rna|ch3|cnv|protein] (default: all)
    --skip_*              Flags to skip specific integration steps
    --emit_qc             Per-case QC to out_dir/qc/{case-id}_qc.tsv
    --keep-join-cols      [off|on] On keeps the linker mapping in the final output csv (default = off)

"""
import pandas as pd

logging.basicConfig(level=logging.INFO, format='[multiomic_integration] %(levelname)s %(message)s')

# ----------------------------- Utilities -----------------------------

def _ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)

def _safe_to_numeric(s: pd.Series) -> pd.Series:
    return pd.to_numeric(s, errors="coerce")

def _norm_ensg_ser(s: pd.Series) -> pd.Series:
    """Extract versionless Ensembl IDs (ENSG#########) from strings."""
    return s.astype(str).str.extract(r'(ENSG\d+)', expand=False)

def _read_any_table(path: Path, sep: Optional[str] = None, names: Optional[List[str]] = None) -> pd.DataFrame:
    """Robust table reader for CSV/TSV/TXT; falls back to python engine on failure."""
    if sep is None:
        sep = "," if path.suffix.lower() == ".csv" else "\t"
    try:
        return pd.read_csv(path, sep=sep, engine="c", dtype=str, low_memory=False,
                           header=None if names else "infer", names=names)
    except Exception:
        return pd.read_csv(path, sep=sep, engine="python", dtype=str,
                           header=None if names else "infer", names=names)


def _normalize_case_sample_columns(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    lc = {c.lower().strip(): c for c in df.columns}
    case_c = lc.get("case id") or lc.get("case") or lc.get("case submitter id")
    samp_c = lc.get("sample id") or lc.get("sample") or lc.get("sample submitter id")
    # Explode comma-separated lists
    if case_c and df[case_c].astype(str).str.contains(",").any():
        df[case_c] = df[case_c].astype(str).str.split(",").map(lambda xs: [x.strip() for x in xs])
        df = df.explode(case_c, ignore_index=True)
        df[case_c] = df[case_c].astype(str).str.strip()
    if samp_c and df[samp_c].astype(str).str.contains(",").any():
        df[samp_c] = df[samp_c].astype(str).str.split(",").map(lambda xs: [x.strip() for x in xs])
        df = df.explode(samp_c, ignore_index=True)
        df[samp_c] = df[samp_c].astype(str).str.strip()
    return df


def _read_first_two_col_table(path: Path, gz_ok: bool = True) -> pd.DataFrame:
    """Read a two-column TSV; supports .gz if gz_ok."""
    if path.suffix == ".gz" and gz_ok:
        with gzip.open(path, "rt") as fh:
            return pd.read_csv(fh, sep="\t", header=None, names=["col1", "col2"], dtype=str)
    else:
        return pd.read_csv(path, sep="\t", header=None, names=["col1", "col2"], dtype=str)

def _try_case_file_in_dir(mod_dir: Path, case_id: str) -> Optional[Path]:
    """Return {case}.{csv|tsv|txt} under mod_dir if present."""
    for ext in (".csv", ".tsv", ".txt"):
        p = mod_dir / f"{case_id}{ext}"
        if p.exists():
            return p
    return None

def _find_in_subdir_by_patterns(root: Path, patterns: List[str]) -> Optional[Path]:
    """Find first file beneath root whose name contains ALL lowercase substrings in patterns."""
    if not root or not root.exists():
        return None
    for p in root.rglob("*"):
        if p.is_file():
            name = p.name.lower()
            if all(pat in name for pat in patterns):
                return p
    return None

def _find_in_subdir_by_patterns(root: Path, pattern_groups: List[Tuple[str, ...]]) -> Optional[Path]:
    """
    Return the first file under `root` whose filename contains ALL tokens in ANY group.
    Each group is ANDed; groups are ORed.
    Tokens inside a group may be strings, or iterables of strings (treated as OR within that position).
    """
    if not root or not root.exists():
        return None

    for p in root.rglob("*"):
        if not p.is_file():
            continue
        name = p.name.lower()

        def token_ok(tok) -> bool:
            # Accept nested tuples/lists/sets: any of them is OK for this token position
            if isinstance(tok, (tuple, list, set)):
                return any(isinstance(t, str) and t in name for t in tok)
            return isinstance(tok, str) and tok in name

        for group in pattern_groups:
            # Accept accidental bare-string groups
            if isinstance(group, str):
                group = (group,)
            elif not isinstance(group, (tuple, list, set)):
                group = (str(group),)

            if all(token_ok(tok) for tok in group):
                return p
    return None

def _resolve_normal_dir(normal_root: Path, modality: str) -> Path:
    """
    Accept either a parent (contains 'normal/') or the 'normal' directory itself.
    Also tolerate 'cn' vs 'cnv'.
    """
    candidates = []
    if normal_root.name.lower() == "normal":
        candidates += [normal_root / modality]
        if modality == "cn":
            candidates += [normal_root / "cnv"]
    candidates += [normal_root / "normal" / modality, normal_root / modality]
    if modality == "cn":
        candidates += [normal_root / "normal" / "cnv", normal_root / "cnv"]

    for c in candidates:
        if c.exists():
            return c
    return candidates[0]

# ------------------------ Reference resolution -----------------------

def _pick_best_subdir(tmp_root: Path) -> Path:
    subs = [d for d in tmp_root.iterdir() if d.is_dir()]
    for d in subs:
        if d.name.lower() == "reference":
            return d
    return subs[0] if subs else tmp_root

def resolve_reference_dir(ref_dir: Optional[str], ref_zip: Optional[str]) -> Optional[Path]:
    """Return extracted reference directory from --ref_dir or --ref_zip (zip is extracted to temp)."""
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

def _find_ref_file(ref_dir: Optional[Path], *name_variants: str) -> Optional[Path]:
    """
    Find a reference file by trying several name variants.
    Each variant may be a stem (no extension) or a full filename.
    """
    if ref_dir is None:
        return None
    for name in name_variants:
        direct = ref_dir / name
        if direct.exists():
            return direct
    # search whole tree for any matching stem or name prefix
    targets = set(name_variants)
    stems = set(Path(n).stem for n in name_variants)
    for p in ref_dir.rglob("*"):
        if p.is_file():
            if p.name in targets or p.stem in stems:
                return p
    return None

# --------------------------- Case manifest ---------------------------

def load_case_ids(manifest_path: str) -> List[str]:
    """
    Load case IDs from a manifest: .txt (one per line, '#' comments ok) or CSV/TSV (first column).
    """
    p = Path(manifest_path)
    if not p.exists():
        raise FileNotFoundError(f"Manifest not found: {p}")

    def dedup_keep_order(seq):
        seen, out = set(), []
        for x in seq:
            if x not in seen:
                seen.add(x); out.append(x)
        return out

    if p.suffix.lower() == ".txt":
        ids = []
        for ln in p.read_text(encoding="utf-8", errors="replace").splitlines():
            ln = ln.strip().lstrip("\ufeff")
            if not ln or ln.startswith("#"):
                continue
            ids.append(ln)
        ids = dedup_keep_order(ids)
        if not ids:
            raise ValueError(f"No case IDs found in {p}")
        return ids

    sepx = "\t" if p.suffix.lower() == ".tsv" else ","
    for hdr in (None, 0):
        try:
            s = pd.read_csv(p, sep=sepx, dtype=str, usecols=[0], header=hdr).iloc[:, 0]
            s = s.dropna().astype(str).str.strip().str.lstrip("\ufeff")
            header_like = {"case-id","case_id","case id","sample-id","sample_id","sample id","id"}
            s = s[~s.str.lower().isin(header_like)]
            ids = [x for x in s.tolist() if x]
            ids = dedup_keep_order(ids)
            if ids:
                return ids
        except Exception:
            continue
    raise ValueError(f"Could not parse case IDs from {p}")

# --------- Auto-detect per-modality manifests & path resolution ------

def _normalize_col_map(df: pd.DataFrame) -> dict:
    import re as _re
    return {_re.sub(r"\s+", " ", c.lower().replace("-", " ").replace("_", " ").strip()): c for c in df.columns}

def _autodetect_manifest(mod_dir: Path, modality: str) -> Optional[Path]:
    """
    If user didn't pass a manifest path, pick a top-level file in mod_dir that
    looks like a manifest (has 'Case ID' and one of 'File ID'/'File Name'/'Path').
    Prefer filenames that include the modality (e.g., 'rna','ch3','cn','cnv') and 'tumor'/'normal'.
    """
    if not mod_dir.exists():
        return None
    top_files = [p for p in mod_dir.iterdir() if p.is_file() and p.suffix.lower() in {".tsv",".txt",".csv"}]
    if not top_files:
        return None

    def score_name(p: Path) -> int:
        name = p.name.lower()
        score = 0
        if modality == "rna":
            for tok in ("rna","htseq","counts"): 
                if tok in name: score += 2
        elif modality == "ch3":
            for tok in ("ch3","methyl","beta","illumina"): 
                if tok in name: score += 2
        elif modality == "cn":
            for tok in ("cn","cnv","copy","gene_level"): 
                if tok in name: score += 2
        for tok in ("tumor","normal"):
            if tok in name: score += 1
        return score

    candidates: List[Tuple[int, Path]] = []
    for f in top_files:
        try:
            df = _read_any_table(f, sep=None)
            norm = _normalize_col_map(df)
            has_case = any(k in norm for k in ("case id","case","sample id","sample"))
            has_file = any(k in norm for k in ("file id","id","file name","filename","name","path","file path","filepath","relpath","full path"))
            if has_case and has_file:
                candidates.append((score_name(f), f))
        except Exception:
            continue
    if not candidates:
        return None
    candidates.sort(key=lambda t: t[0], reverse=True)
    best = candidates[0][1]
    logging.info(f"[AUTO] Using manifest {best} for modality={modality}")
    return best

def _load_manifest_df(mod_dir: Path, manifest: Optional[str], default_name: str, modality: str) -> Optional[pd.DataFrame]:
    """
    Load a case→file manifest for a modality.
    Accepts absolute path or path relative to mod_dir, or auto-detects a top-level file in mod_dir.
    Normalizes to columns: 'Case ID' + any of ['File ID', 'File Name', 'Path'].
    """
    candidates: List[Path] = []
    if manifest:
        mp = Path(manifest)
        candidates += [mp, mod_dir / mp.name, mod_dir / str(manifest)]
    candidates += [mod_dir / default_name]

    man_path = next((p for p in candidates if p and p.exists()), None)
    if man_path is None:
        man_path = _autodetect_manifest(mod_dir, modality)
        if man_path is None:
            return None

    df = _read_any_table(man_path, sep="\t" if man_path.suffix.lower() in {".tsv",".txt"} else None)
    df = _normalize_case_sample_columns(df)
    norm = _normalize_col_map(df)

    def pick(*keys):
        for k in keys:
            if k in norm:
                return norm[k]
        return None

    case_c = pick("case id", "case", "sample id", "sample")
    fileid_c = pick("file id", "id")
    fname_c = pick("file name", "filename", "name")
    path_c  = pick("path", "file path", "filepath", "relpath", "full path")

    if not case_c or not (fileid_c or fname_c or path_c):
        logging.warning(f"Manifest {man_path} missing required columns (need Case ID and one of File ID / File Name / Path).")
        return None

    out = pd.DataFrame({"Case ID": df[case_c].astype(str).str.strip()})
    if fileid_c: out["File ID"] = df[fileid_c].astype(str).str.strip()
    if fname_c:  out["File Name"] = df[fname_c].astype(str).str.strip()
    if path_c:   out["Path"] = df[path_c].astype(str).str.strip()
    return out.dropna(subset=["Case ID"])

def _resolve_data_root(mod_dir: Path, row: pd.Series) -> Optional[Path]:
    """
    Given a manifest row, return either:
      - a directory to search (e.g., <mod_dir>/<File ID>/),
      - a direct file path if 'Path' or 'File Name' point to a file,
      - or best-effort match by basename under mod_dir.
    """
    def _norm_join(base: Path, p: str) -> Path:
        q = Path(str(p).strip())
        if not q.is_absolute():
            q = (base / q).resolve()
        return q

    # 1) Explicit Path
    path_val = row.get("Path")
    if pd.notna(path_val):
        p = _norm_join(mod_dir, str(path_val))
        if p.exists(): 
            return p
        if p.name:
            hit = next((h for h in mod_dir.rglob(p.name)), None)
            if hit: 
                return hit
        try:
            alt = mod_dir.joinpath(*Path(str(path_val)).parts)
            if alt.exists():
                return alt
        except Exception:
            pass

    # 2) GDC-style File ID directory
    fid = row.get("File ID")
    if pd.notna(fid):
        candidate = mod_dir / str(fid).strip()
        if candidate.exists():
            return candidate

    # 3) File Name (may include subdirs)
    fname = row.get("File Name")
    if pd.notna(fname):
        fname = str(fname).strip()
        if "/" in fname or "\\" in fname:
            p2 = _norm_join(mod_dir, fname)
            if p2.exists(): 
                return p2
        hit = next((h for h in mod_dir.rglob(Path(fname).name)), None)
        if hit: 
            return hit
    return None

def _emit_qc_counts(case_id: str, base_df):
    # Tumor counts
    rna_t = int(base_df["RNA_Count"].notna().sum()) if "RNA_Count" in base_df.columns else 0
    cnv_t = (int(base_df["CNV_Count"].notna().sum()) if "CNV_Count" in base_df.columns else          (int(base_df["copy_number"].notna().sum()) if "copy_number" in base_df.columns else 0))
    pro_t = int(base_df["SEQ"].notna().sum() or base_df["Protein_SEQ"].notna().sum()) if ("SEQ" in base_df.columns or "Protein_SEQ" in base_df.columns) else 0
    ch3_t = int(base_df["CH3_Beta"].notna().sum()) if "CH3_Beta" in base_df.columns else 0

    if rna_t: logging.info(f"[RNA] {case_id}: matched RNA_Count for {rna_t} rows")
    if cnv_t: logging.info(f"[CNV] {case_id}: matched CNV_Count for {cnv_t} rows")
    if pro_t: logging.info(f"[PROTEIN] {case_id}: matched NP for {pro_t} rows")
    if ch3_t: logging.info(f"[CH3] {case_id}: matched CH3_Beta for {ch3_t} rows")

    # Normal counts
    rna_n = int(base_df["Normal_RNA_Count"].notna().sum()) if "Normal_RNA_Count" in base_df.columns else 0
    cnv_n = (int(base_df["Normal_CNV_Count"].notna().sum()) if "Normal_CNV_Count" in base_df.columns else          (int(base_df["Normal_copy_number"].notna().sum()) if "Normal_copy_number" in base_df.columns else 0))
    pro_n = int(base_df["Normal_SEQ"].notna().sum()) if "Normal_SEQ" in base_df.columns else 0
    ch3_n = int(base_df["Normal_CH3_Beta"].notna().sum()) if "Normal_CH3_Beta" in base_df.columns else 0

    if rna_n: logging.info(f"[RNA:NORMAL] {case_id}: matched Normal_RNA_Count for {rna_n} rows")
    if cnv_n: logging.info(f"[CNV:NORMAL] {case_id}: matched Normal_CNV_Count for {cnv_n} rows")
    if pro_n: logging.info(f"[PROTEIN:NORMAL] {case_id}: matched Normal_SEQ for {pro_n} rows")
    if ch3_n: logging.info(f"[CH3:NORMAL] {case_id}: matched Normal_CH3_Beta for {ch3_n} rows")


# ------------------------------ DNA ---------------------------------

def _add_ENSG_columns(df: pd.DataFrame, case_id: str) -> pd.DataFrame:
    """Ensure ENSGene & ENSGene_core exist; pull ENSG from INFO if available."""
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
    # Always add versionless
    if "ENSGene" in df.columns:
        df["ENSGene_core"] = _norm_ensg_ser(df["ENSGene"])
    elif "INFO" in df.columns:
        df["ENSGene_core"] = _norm_ensg_ser(df["INFO"])
    else:
        df["ENSGene_core"] = pd.Series(index=df.index, dtype=object)
    return df

def load_and_enhance_dna(case_id: str, dna_dir: Path, out_dir: Path,
                         ref_dir: Optional[Path], dna_rows: str) -> Path:
    """Read per-case DNA, add ENSG columns, add maps, write out_dir/dna/{case}.csv."""
    inp = _try_case_file_in_dir(dna_dir, case_id)
    if inp is None:
        tried = [str(dna_dir / f"{case_id}{ext}") for ext in (".csv", ".tsv", ".txt")]
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
    out_dna = out_dir / "dna"
    _ensure_dir(out_dna)
    out_fp = out_dna / f"{case_id}.csv"
    df.to_csv(out_fp, index=False)
    return out_fp

# --------------------------- Reference maps --------------------------

def _sym_file(ref_dir: Optional[Path]) -> Optional[Path]:
    """Find symbol↔ENSG map: esng_gene-sym(.txt) or ensg_gene-sym(.txt)."""
    if ref_dir is None:
        return None
    return (_find_ref_file(ref_dir, "esng_gene-sym.txt", "esng_gene-sym")
            or _find_ref_file(ref_dir, "ensg_gene-sym.txt", "ensg_gene-sym"))

def integrate_gene_annotations(base_df: pd.DataFrame, ref_dir: Optional[Path], case_id: str) -> pd.DataFrame:
    """Attach Gene symbol and name via reference (esng/ensg_gene-sym, gene-sym_name)."""
    if ref_dir is None:
        return base_df
    df = base_df.copy()

    sym_fp  = _sym_file(ref_dir)
    name_fp = _find_ref_file(ref_dir, "gene-sym_name.txt", "gene-sym_name")

    if sym_fp:
        m = _read_any_table(sym_fp, sep="\t", names=["Gene", "Ensembl"]).drop_duplicates()
        m["Gene"] = m["Gene"].astype(str).str.strip()
        m["Ensembl_core"] = _norm_ensg_ser(m["Ensembl"])
        m = m.drop_duplicates(subset=["Ensembl_core"], keep="first")
        df = df.merge(m[["Gene", "Ensembl_core"]], how="left",
                      left_on="ENSGene_core", right_on="Ensembl_core")
    else:
        logging.warning(f"[DNA] {case_id}: Gene symbol map not found in reference")

    if name_fp and "Gene" in df.columns:
        n = _read_any_table(name_fp, sep="\t", names=["symbol", "name"]).drop_duplicates()
        n["symbol"] = n["symbol"].astype(str).str.strip()
        n["name"] = n["name"].astype(str).str.strip()
        n = n.drop_duplicates(subset=["symbol"], keep="first")
        df = df.merge(n, how="left", left_on="Gene", right_on="symbol")

    # Cleanup temp columns from this step
    for col in ("Ensembl_core", "symbol"):
        if col in df.columns:
            df = df.drop(columns=[col])
    if "name" in df.columns and "Gene_Name" not in df.columns:
        df = df.rename(columns={"name": "Gene_Name"})
    return df

def integrate_uniprot_np(base_df: pd.DataFrame, ref_dir: Optional[Path], case_id: str) -> pd.DataFrame:
    """Attach UniProt and NP/XP (if available)."""
    if ref_dir is None:
        return base_df
    df = base_df.copy()
    ue_fp  = _find_ref_file(ref_dir, "uni-ensg_all.txt")
    unp_fp = _find_ref_file(ref_dir, "uni-np_all.txt")
    if not ue_fp:
        logging.warning(f"[DNA] {case_id}: UniProt→ENSG map not found")
        return df

    ue = _read_any_table(ue_fp, sep="\t", names=["UniProt", "ENSG"]).drop_duplicates()
    ue["To_core"] = _norm_ensg_ser(ue["ENSG"])
    ens_to_up = (ue.groupby("To_core")["UniProt"]
                   .apply(lambda s: sorted(set(str(x).strip() for x in s if pd.notna(x))))
                   .reset_index(name="UniProt_list"))
    df = df.merge(ens_to_up, how="left", left_on="ENSGene_core", right_on="To_core").drop(columns=["To_core"])

    if unp_fp and "UniProt_list" in df.columns:
        unp = _read_any_table(unp_fp, sep="\t", names=["UniProt", "NP_ID"]).drop_duplicates()
        unp["UniProt"] = unp["UniProt"].astype(str).str.strip()
        up_to_nps = (unp.groupby("UniProt")["NP_ID"]
                        .apply(lambda s: sorted(set(str(x).strip() for x in s if pd.notna(x))))
                        .to_dict())

        def _collect_nps(up_list):
            if not isinstance(up_list, list):
                return None
            bag = []
            for up in up_list:
                bag.extend(up_to_nps.get(up, []))
            bag = sorted(set(bag))
            return bag if bag else None

        df["NP_list"] = df["UniProt_list"].map(_collect_nps)
        df["UniProt"] = df["UniProt_list"].map(lambda xs: ";".join(xs) if isinstance(xs, list) and xs else None)
        df["NP"] = df["NP_list"].map(lambda xs: ";".join(xs) if isinstance(xs, list) and xs else None)
        df = df.drop(columns=["UniProt_list", "NP_list"])
    elif "UniProt_list" in df.columns:
        df["UniProt"] = df["UniProt_list"].map(lambda xs: ";".join(xs) if isinstance(xs, list) and xs else None)
        df = df.drop(columns=["UniProt_list"])
    return df

# ------------------------------- RNA --------------------------------

def _find_rna_counts_file(root_or_file: Path) -> Optional[Path]:
    """
    Prefer common gene counts outputs; otherwise fallback to any table that looks like ENSG counts.
    """
    if root_or_file.is_file():
        return root_or_file

    # OR of AND-groups (tokens are ANDed inside each tuple)
    patterns = [
        ("gene", "quantification"),        # your header example
        ("rsem", "genes"),
        ("htseq", "counts"),
        ("readspergene",),                 # STAR ReadsPerGene.out.tab
        ("counts",),
        (".tsv",), (".txt",), (".tsv.gz",), (".txt.gz",),
    ]
    hit = _find_in_subdir_by_patterns(root_or_file, patterns)
    if hit:
        return hit

    # Fallback: scan small text/tsv[.gz] files and pick first that looks like ENSG table
    for p in root_or_file.rglob("*"):
        if not p.is_file():
            continue
        if p.suffix.lower() in {".tsv", ".txt"} or (p.suffix.lower() == ".gz" and p.name.endswith((".tsv.gz", ".txt.gz"))):
            try:
                opener = gzip.open if p.suffix == ".gz" else open
                with opener(p, "rt") as fh:
                    head = pd.read_csv(fh, sep="\t", nrows=50, dtype=str, comment="#")
                cols = [c.lower() for c in head.columns]
                if any(c.startswith("gene") or "ensembl" in c for c in cols):
                    return p
            except Exception:
                continue
    return None

def _rna_from_manifest(rna_dir: Path, case_id: str, rna_manifest: Optional[str]) -> Optional[pd.DataFrame]:
    man = _load_manifest_df(rna_dir, rna_manifest, "gdc-rna.tsv", modality="rna")
    if man is None:
        return None
    rows = man[man["Case ID"] == case_id]
    if rows.empty:
        logging.warning(f"[CNV] {case_id}: no manifest row matched; check Case ID exact text and commas.")
        return None

    root_or_file = _resolve_data_root(rna_dir, rows.iloc[0])
    if root_or_file is None:
        logging.warning(f"[RNA] {case_id}: cannot resolve path from manifest.")
        return None

    f = _find_rna_counts_file(root_or_file)
    if f is None or not f.exists():
        logging.warning(f"[RNA] {case_id}: no RNA counts file found under {root_or_file}")
        return None

    logging.info(f"[RNA] {case_id}: using {f.name}")

    # Try multi-column gene-quantification first (your format)
    try:
        df0 = pd.read_csv(f, sep="\t", dtype=str, comment="#")
        lc = {c.lower(): c for c in df0.columns}

        gene_c = lc.get("gene_id") or lc.get("ensembl") or lc.get("ensembl_gene_id") or lc.get("gene")
        # Prefer raw counts; fall back to TPM/FPKM if needed
        count_c = (lc.get("unstranded") or lc.get("raw_counts") or lc.get("counts") or
                   lc.get("stranded_first") or lc.get("stranded_second") or
                   lc.get("tpm_unstranded") or lc.get("tpm") or
                   lc.get("fpkm_uq_unstranded") or lc.get("fpkm"))

        if gene_c and count_c:
            df = df0[[gene_c, count_c]].rename(columns={gene_c: "Ensembl_full", count_c: "Count"}).copy()
            # Keep only real genes (drop N_unmapped, etc.)
            df = df[df["Ensembl_full"].astype(str).str.match(r"^ENSG", na=False)]
            df["Ensembl_core"] = _norm_ensg_ser(df["Ensembl_full"])  # strips .15
            df["Count"] = _safe_to_numeric(df["Count"])
            df = df.drop_duplicates(subset=["Ensembl_full"], keep="first")
            return df
    except Exception as e:
        logging.debug(f"[RNA] {case_id}: multi-column parse failed: {e}")

    # Fallback: HTSeq-style 2-column table
    df2 = _read_first_two_col_table(f, gz_ok=True).rename(columns={"col1": "Ensembl_full", "col2": "Count"})
    df2 = df2[df2["Ensembl_full"].astype(str).str.match(r"^ENSG", na=False)]
    df2["Ensembl_core"] = _norm_ensg_ser(df2["Ensembl_full"])
    df2["Count"] = _safe_to_numeric(df2["Count"])
    df2 = df2.drop_duplicates(subset=["Ensembl_full"], keep="first")
    return df2


def integrate_rna(base_df: pd.DataFrame, rna_dir: Path, case_id: str,
                  ensg_join_mode: str, synth_overlap: bool, rna_manifest: Optional[str]) -> pd.DataFrame:
    df = base_df.copy()
    rna = _rna_from_manifest(rna_dir, case_id, rna_manifest)
    if rna is None or rna.empty:
        return df
    left_key = "ENSGene_core" if ensg_join_mode == "core" else "ENSGene"
    right_key = "Ensembl_core" if ensg_join_mode == "core" else "Ensembl_full"
    if synth_overlap:
        need = set(df[left_key].dropna().unique()) - set(rna[right_key].dropna().unique())
        if need:
            add = pd.DataFrame({
                "Ensembl_full": [x if ensg_join_mode == "exact" else
                                 (f"{x}.1" if not re.search(r'\.\d+$', str(x)) else str(x)) for x in need],
                "Ensembl_core": list(need), "Count": 100
            })
            rna = pd.concat([rna, add], ignore_index=True)
    df = df.merge(rna[[right_key, "Count"]], how="left", left_on=left_key, right_on=right_key)
    if right_key in df.columns:
        df = df.drop(columns=[right_key])
    logging.info(f"[RNA] {case_id}: matched Count for {int(df['Count'].notna().sum())} rows (join={ensg_join_mode})")
    return df

# ------------------------------- CNV --------------------------------

def _find_cnv_gene_level_file(root_or_file: Path) -> Optional[Path]:
    """Looser discovery for CNV (gene-level)."""
    if root_or_file.is_file():
        return root_or_file

    # Filename heuristics (OR of AND-groups)
    patterns = [
        ("gene", "copy", "number"),
        ("gene", "level", "copy"),
        ("gene_level", "copy_number"), ("gene_level.copy_number_variation",),
        ("copy_number_variation",),
        ("gene", "cnv"),
        ("copy_ratio",),
        ("segment_mean",),
        ("log2", "ratio"),
    ]
    hit = _find_in_subdir_by_patterns(root_or_file, patterns)
    if hit:
        return hit

    # Column-based fallback: find a TSV/CSV that has a gene id and a copy-number-ish column
    for p in list(root_or_file.rglob("*.tsv")) + list(root_or_file.rglob("*.csv")) + list(root_or_file.rglob("*.txt")):
        try:
            head = pd.read_csv(p, sep="\t" if p.suffix.lower() != ".csv" else ",", nrows=30, dtype=str)
            lc = {c.lower().replace("-", " ").replace("_", " ").strip(): c for c in head.columns}
            has_gene = any(k in lc for k in ("gene id","ensembl","ensembl gene id","gene","symbol"))
            has_cn   = any(k in lc for k in ("copy number","cn","cn value","total copy number","integer copy number",
                                             "min copy number","max copy number","copy ratio","segment mean","log2"))
            if has_gene and has_cn:
                return p
        except Exception:
            continue
    return None

def _cn_from_manifest(cn_dir: Path, case_id: str, cn_manifest: Optional[str], ref_dir: Optional[Path]=None) -> Optional[pd.DataFrame]:
    man = _load_manifest_df(cn_dir, cn_manifest, "gdc-cn.tsv", modality="cn")
    if man is None:
        return None

    rows = man[man["Case ID"] == case_id]
    if rows.empty:
        logging.warning(f"[CNV] {case_id}: no manifest row matched; check Case ID exact text and commas.")
        return None

    root_or_file = _resolve_data_root(cn_dir, rows.iloc[0])
    if root_or_file is None:
        logging.warning(f"[CNV] {case_id}: cannot resolve path from manifest.")
        return None

    f = _find_cnv_gene_level_file(root_or_file)
    if f is None or not f.exists():
        logging.warning(f"[CNV] {case_id}: no gene-level CNV file found under {root_or_file}")
        return None

    sep = "," if f.suffix.lower() == ".csv" else "\t"
    logging.info(f"[CNV] {case_id}: using {f.name}")
    df0 = pd.read_csv(f, sep=sep, dtype=str, comment="#")
    # Fast-path for common CPTAC-style gene-level CNV table
    # Expected columns (example): gene_id, gene_name, chromosome, start, end, copy_number, min_copy_number, max_copy_number
    cols_lower = [c.lower().strip() for c in df0.columns]
    if {"gene_id","copy_number"}.issubset(set(cols_lower)):
        # Map to canonical names
        col_map = {c.lower().strip(): c for c in df0.columns}
        gid_col = col_map["gene_id"]
        cn_col  = col_map["copy_number"]
        # Build standardized CN table
        cn_out = pd.DataFrame({
            "Ensembl_full": df0[gid_col].astype(str),
            "Ensembl_core": _norm_ensg_ser(df0[gid_col]),
            "copy_number": _safe_to_numeric(df0[cn_col]),
        })
        cn_out = cn_out.dropna(subset=["Ensembl_core"]).drop_duplicates(subset=["Ensembl_core"], keep="first")
        return cn_out

    lc = {c.lower().replace("-", " ").replace("_", " ").strip(): c for c in df0.columns}

    # Identify gene column
    gid_c = lc.get("gene id") or lc.get("ensembl") or lc.get("ensembl gene id")
    symbol_c = lc.get("gene") or lc.get("symbol")

    # Identify copy-number-like column
    cn_c  = (lc.get("copy number") or lc.get("cn") or lc.get("cn value") or
             lc.get("total copy number") or lc.get("integer copy number") or
             lc.get("copy ratio") or lc.get("segment mean") or lc.get("log2"))

    if not (gid_c or symbol_c) or not cn_c:
        logging.warning(f"[CNV] {case_id}: expected gene_id/symbol and copy_number-like column in {f.name}")
        return None

    # Normalize CN values; if segment/log2 provided, approximate CN = 2 * 2^(log2ratio)
    cn_series = pd.to_numeric(df0[cn_c], errors="coerce")
    if cn_c in (lc.get("copy ratio"), lc.get("segment mean"), lc.get("log2")):
        # Assume these represent log2 ratio; convert
        cn_series = 2.0 * (2.0 ** cn_series)

    out = pd.DataFrame({"copy_number": cn_series})

    if gid_c and df0[gid_c].astype(str).str.contains(r"ENSG", na=False).any():
        out["Ensembl_full"] = df0[gid_c].astype(str)
        out["Ensembl_core"] = _norm_ensg_ser(out["Ensembl_full"])
    else:
        # Use gene symbols → ENSG map
        if symbol_c is None:
            logging.warning(f"[CNV] {case_id}: no ENSG and no symbol column found in {f.name}")
            return None
        sym_fp = _sym_file(ref_dir) if ref_dir else None
        if sym_fp is None:
            logging.warning(f"[CNV] {case_id}: symbol-based CNV but no reference symbol map available")
            return None
        sym = _read_any_table(sym_fp, sep="\t", names=["Gene", "Ensembl"]).drop_duplicates()
        sym["Gene"] = sym["Gene"].astype(str).str.strip()
        sym["Ensembl_core"] = _norm_ensg_ser(sym["Ensembl"])
        sym = sym.drop_duplicates(subset=["Gene"], keep="first")
        tmp = df0[[symbol_c]].rename(columns={symbol_c:"Gene"}).copy()
        out = pd.concat([tmp, out], axis=1).merge(sym[["Gene","Ensembl_core"]], how="left", on="Gene")
        out = out.drop(columns=["Gene"])

    out = out.dropna(subset=["Ensembl_core"])
    out["copy_number"] = _safe_to_numeric(out["copy_number"])
    out = out.drop_duplicates(subset=["Ensembl_core"], keep="first")
    return out

def integrate_cnv(base_df: pd.DataFrame, cn_dir: Path, case_id: str,
                  ensg_join_mode: str, synth_overlap: bool, cn_manifest: Optional[str],
                  ref_dir: Optional[Path]=None) -> pd.DataFrame:
    df = base_df.copy()
    cn = _cn_from_manifest(cn_dir, case_id, cn_manifest, ref_dir=ref_dir)
    if cn is None or cn.empty:
        return df
    left_key = "ENSGene_core" if ensg_join_mode == "core" else "ENSGene"
    right_key = "Ensembl_core" if ensg_join_mode == "core" else "Ensembl_full"
    if synth_overlap:
        need = set(df[left_key].dropna().unique()) - set(cn[right_key].dropna().unique())
        if need:
            add = pd.DataFrame({
                "Ensembl_full": [x if ensg_join_mode == "exact" else
                                 (f"{x}.1" if not re.search(r'\.\d+$', str(x)) else str(x)) for x in need],
                "Ensembl_core": list(need), "copy_number": 2
            })
            cn = pd.concat([cn, add], ignore_index=True)
    lk = left_key; rk = right_key
    luniq = int(df[lk].dropna().nunique())
    runiq = int(cn[rk].dropna().nunique())
    df = df.merge(cn[[rk, "copy_number"]], how="left", left_on=lk, right_on=rk)
    matched = int(df["copy_number"].notna().sum())
    logging.info(f"[CNV] {case_id}: join keys left={luniq} vs right={runiq}; matched copy_number in {matched} rows")
    logging.info(f"[CNV] {case_id}: matched copy_number for {matched} rows")
    if right_key in df.columns:
        df = df.drop(columns=[right_key])
    logging.info(f"[CNV] {case_id}: matched copy_number for {int(df['copy_number'].notna().sum())} rows (join={ensg_join_mode})")
    return df

# ------------------------------- CH3 --------------------------------

def _ch3_from_manifest(ch3_dir: Path, case_id: str, ch3_manifest: Optional[str]) -> Optional[pd.DataFrame]:
    man = _load_manifest_df(ch3_dir, ch3_manifest, "gdc-ch3.tsv", modality="ch3")
    if man is None:
        return None
    rows = man[man["Case ID"] == case_id]
    if rows.empty:
        logging.warning(f"[CNV] {case_id}: no manifest row matched; check Case ID exact text and commas.")
        return None
    root_or_file = _resolve_data_root(ch3_dir, rows.iloc[0])
    if root_or_file is None:
        logging.warning(f"[CH3] {case_id}: cannot resolve path from manifest.")
        return None

    candidate = root_or_file if root_or_file.is_file() else None
    if candidate is None:
        # try raw level3betas (no header)
        for p in root_or_file.rglob("*.txt"):
            try:
                head = pd.read_csv(p, sep="\t", dtype=str, nrows=1, header=None)
                if head.shape[1] >= 2:
                    df = pd.read_csv(p, sep="\t", dtype=str, header=None, names=["probe", "beta"])
                    if df.shape[1] >= 2 and df["probe"].str.startswith("cg").any():
                        return df.assign(beta=_safe_to_numeric(df["beta"]))
            except Exception:
                continue
        # fallback: headered TSVs
        for p in root_or_file.rglob("*.tsv"):
            try:
                head = pd.read_csv(p, sep="\t", dtype=str, nrows=16)
                lc = {c.lower(): c for c in head.columns}
                if (lc.get("cg") or lc.get("ilmnid") or lc.get("composite element ref") or lc.get("probe")) and \
                   (lc.get("beta") or lc.get("beta_value") or lc.get("beta value")):
                    candidate = p; break
            except Exception:
                continue
    if candidate is None:
        logging.warning(f"[CH3] {case_id}: could not locate cg/beta table under {root_or_file}")
        return None

    full = pd.read_csv(candidate, sep="\t", dtype=str)
    lc = {c.lower(): c for c in full.columns}
    probe_col = lc.get("probe") or lc.get("cg") or lc.get("ilmnid") or lc.get("composite element ref")
    beta_col  = lc.get("beta") or lc.get("beta_value") or lc.get("beta value")
    if probe_col is None or beta_col is None:
        logging.warning(f"[CH3] {case_id}: unable to identify probe/beta columns in {candidate.name}; "
                        f"found probe_col={probe_col}, beta_col={beta_col}. "
                        f"Available columns (first 10): {list(full.columns[:10])}")
        return None
    out = full[[probe_col, beta_col]].rename(columns={probe_col: "probe", beta_col: "beta"})
    out["beta"] = _safe_to_numeric(out["beta"])
    return out.dropna(subset=["probe"])

def _load_ch3_map(ch3_map: Optional[str], ref_dir: Optional[Path],
                  ch3_probe_col: Optional[str], ch3_ensg_col: Optional[str], ch3_symbol_col: Optional[str]) -> Optional[pd.DataFrame]:
    """Load a probe→ENSG map via path or by searching the reference for common filenames."""
    path: Optional[Path] = None
    if ch3_map:
        path = Path(ch3_map)
    else:
        for name in ["ch3.csv", "ch3_map.csv", "probe-ensg.tsv", "probe_ensg.tsv",
                     "illumina_probes_ensg.tsv", "ch3_probe_map.tsv", "ch3_map.tsv"]:
            guess = _find_ref_file(ref_dir, name) if ref_dir else None
            if guess and guess.exists():
                path = guess; break
    if path is None or not path.exists():
        return None

    sep = "," if path.suffix.lower() == ".csv" else "\t"
    df = _read_any_table(path, sep=sep)
    lc = {c.lower(): c for c in df.columns}

    probe_c = (ch3_probe_col
               or lc.get("ilmnid")          # IlmnID
               or lc.get("probe")
               or lc.get("cg")
               or lc.get("composite element ref")
               or lc.get("illmnid"))        # tolerate misspelling

    ensg_c  = (ch3_ensg_col
               or lc.get("ensg")
               or lc.get("ensembl")
               or lc.get("ensembl_gene_id")
               or lc.get("gene_id"))

    sym_c   = (ch3_symbol_col
               or lc.get("ucsc_refgene_name")
               or lc.get("symbol")
               or lc.get("gene")
               or lc.get("gene_symbol"))

    if not probe_c or not (ensg_c or sym_c):
        logging.warning(f"[CH3] Map {path} missing probe+gene columns; "
                        f"provide --ch3_probe_col and one of --ch3_ensg_col/--ch3_symbol_col.")
        return None

    m = df.rename(columns={probe_c: "probe"}).copy()
    if ensg_c:
        m["ENSGene_core"] = _norm_ensg_ser(m[ensg_c])
    else:
        sym_fp = _sym_file(ref_dir) if ref_dir else None
        if sym_fp is None:
            logging.warning("[CH3] No symbol map (esng/ensg_gene-sym) in reference; cannot derive ENSG from symbol.")
            return None
        sym = _read_any_table(sym_fp, sep="\t", names=["Gene", "Ensembl"]).drop_duplicates()
        sym["Gene"] = sym["Gene"].astype(str).str.strip()
        sym["Ensembl_core"] = _norm_ensg_ser(sym["Ensembl"])
        sym = sym.drop_duplicates(subset=["Gene"], keep="first")
        m = m.rename(columns={sym_c: "Gene"}).merge(
            sym[["Gene", "Ensembl_core"]], how="left", on="Gene"
        ).rename(columns={"Ensembl_core": "ENSGene_core"})

    m = m.dropna(subset=["probe", "ENSGene_core"]).drop_duplicates(subset=["probe"], keep="first")
    return m[["probe", "ENSGene_core"]]

def integrate_ch3(base_df: pd.DataFrame, ch3_dir: Path, case_id: str, ensg_join_mode: str,
                  ch3_map: Optional[str], ref_dir: Optional[Path],
                  ch3_probe_col: Optional[str], ch3_ensg_col: Optional[str], ch3_symbol_col: Optional[str],
                  ch3_agg: str, ch3_manifest: Optional[str]) -> pd.DataFrame:
    df = base_df.copy()
    ch3 = _ch3_from_manifest(ch3_dir, case_id, ch3_manifest)
    if ch3 is None or ch3.empty:
        return df

    m = _load_ch3_map(ch3_map, ref_dir, ch3_probe_col, ch3_ensg_col, ch3_symbol_col)
    if m is None or m.empty:
        logging.warning(f"[CH3] {case_id}: No usable probe→ENSG map; skipping CH3.")
        return df

    ch3m = ch3.merge(m, how="left", on="probe").dropna(subset=["ENSGene_core"])
    if ch3m.empty:
        logging.warning(f"[CH3] {case_id}: probe map didn’t cover this case’s probes; skipping.")
        return df

    # Aggregate beta per gene
    if ch3_agg == "median":
        agg = ch3m.groupby("ENSGene_core", as_index=False)["beta"].median().rename(columns={"beta":"beta_val"})
    elif ch3_agg == "all":
        stats = ch3m.groupby("ENSGene_core")["beta"].agg(['mean','median','min','max','std','count']).reset_index()
        stats = stats.rename(columns={
            'mean':'CH3_Beta_mean','median':'CH3_Beta_median','min':'CH3_Beta_min',
            'max':'CH3_Beta_max','std':'CH3_Beta_std','count':'CH3_Beta_n','ENSGene_core':'ENSGene_core'
        })
        agg = stats
    else:  # mean
        agg = ch3m.groupby("ENSGene_core", as_index=False)["beta"].mean().rename(columns={"beta":"beta_val"})

    # Also collect probe list per gene (helpful for debugging)
    probes = (ch3m.groupby('ENSGene_core')['probe']
              .apply(lambda s: ';'.join(sorted(set(map(str, s)))))
              .reset_index(name='CpG_Probes'))

    left_key = "ENSGene_core" if ensg_join_mode == "core" else "ENSGene"
    df = df.merge(agg, how="left", left_on=left_key, right_on="ENSGene_core")
    if left_key != "ENSGene_core" and "ENSGene_core" in df.columns:
        df = df.drop(columns=["ENSGene_core"])
    df = df.merge(probes, how="left", left_on=left_key, right_on="ENSGene_core")
    if left_key != "ENSGene_core" and "ENSGene_core" in df.columns:
        df = df.drop(columns=["ENSGene_core"])

    if "beta_val" in df.columns:
        logging.info(f"[CH3] {case_id}: matched beta_val for {int(df['beta_val'].notna().sum())} rows (join={ensg_join_mode})")
    else:
        logging.info(f"[CH3] {case_id}: matched aggregated CH3 for {len(agg)} genes (join={ensg_join_mode})")
    return df

# ------------------------------ Protein ------------------------------

def _select_best_protein_rows(prot: pd.DataFrame) -> pd.DataFrame:
    """Keep best row per NP: lowest EV, tie-break by longest SEQ. Drop SEQ length <=3."""
    if "SEQ" not in prot.columns:
        return prot
    prot = prot.copy()
    prot["SEQ_len"] = prot["SEQ"].astype(str).str.len()
    prot = prot[prot["SEQ_len"] > 3]  # drop tiny peptides
    if "EV" in prot.columns:
        prot["EV_num"] = pd.to_numeric(prot["EV"], errors="coerce")
        prot = prot.sort_values(["NP", "EV_num", "SEQ_len"], ascending=[True, True, False])
        prot = prot.drop_duplicates(subset=["NP"], keep="first")
        prot = prot.drop(columns=["EV_num"])
    return prot.drop(columns=["SEQ_len"])

def integrate_protein(base_df: pd.DataFrame, protein_dir: Path, case_id: str, synth_overlap: bool,
                      protein_policy: str = "best") -> pd.DataFrame:
    p = _try_case_file_in_dir(protein_dir, case_id)
    if p is None:
        return base_df
    prot = _read_any_table(p)  # auto sep by extension (.csv/.tsv/.txt)
    prot = prot.drop_duplicates()
    if "NP" not in prot.columns:
        logging.warning(f"[PROTEIN] {case_id}: protein file lacks NP; skipping.")
        return base_df
    if protein_policy == "best":
        prot = _select_best_protein_rows(prot)
    elif protein_policy == "first":
        # Drop SEQ <=3 if present
        if "SEQ" in prot.columns:
            prot = prot[prot["SEQ"].astype(str).str.len() > 3]
        prot = prot.drop_duplicates(subset=["NP"], keep="first")
    else:
        # "all": keep all but drop SEQ <= 3
        if "SEQ" in prot.columns:
            prot = prot[prot["SEQ"].astype(str).str.len() > 3]

    df = base_df.copy()

    # If base has semicolon-joined NP lists, explode them for a precise merge
    if "NP" in df.columns and df["NP"].astype(str).str.contains(";").any():
        df = df.copy()
        df["_orig_idx"] = range(len(df))
        df["NP"] = df["NP"].astype(str).apply(lambda x: [t.strip() for t in x.split(";") if t.strip()] if isinstance(x, str) else [])
        df = df.explode("NP", ignore_index=True)

        merged = df.merge(prot, how="left", on="NP", suffixes=("", "_protein"))
        if protein_policy in ("best", "first"):
            # pick best/first per original row
            if "EV" in merged.columns:
                merged["_ev"] = pd.to_numeric(merged["EV"], errors="coerce")
            else:
                merged["_ev"] = float("inf")
            if "SEQ" in merged.columns:
                merged["_len"] = merged["SEQ"].astype(str).str.len()
            else:
                merged["_len"] = -1
            merged = merged.sort_values(["_orig_idx", "_ev", "_len"], ascending=[True, True, False])
            merged = merged.drop_duplicates(subset=["_orig_idx"], keep="first")
            merged = merged.drop(columns=["_ev","_len","_orig_idx"])
        else:
            merged = merged.drop(columns=["_orig_idx"])
        df = merged
    else:
        df["NP"] = df.get("NP", pd.Series(index=df.index, dtype=object)).astype(str).str.strip()
        df = df.merge(prot, how="left", on="NP", suffixes=("", "_protein"))

    logging.info(f"[PROTEIN] {case_id}: matched NP for {int(df['SEQ'].notna().sum()) if 'SEQ' in df.columns else 0} rows")
    return df

# --------------------------- Dedup & QC ------------------------------

def _final_dedup(df: pd.DataFrame, level: str, key_cols: Optional[List[str]]) -> pd.DataFrame:
    if level == "none":
        return df
    if level == "row":
        return df.drop_duplicates()
    if level == "key" and key_cols:
        cols = [c for c in key_cols if c in df.columns]
        if cols:
            return df.drop_duplicates(subset=cols, keep="first")
    return df

# ---------------------- Normal integration helpers -------------------

def _left_right_keys(ensg_join_mode: str) -> Tuple[str, str]:
    left_key = "ENSGene_core" if ensg_join_mode == "core" else "ENSGene"
    right_key = "Ensembl_core" if ensg_join_mode == "core" else "Ensembl_full"
    return left_key, right_key

def integrate_rna_normal(base_df: pd.DataFrame, normal_root: Path, case_id: str,
                         ensg_join_mode: str, synth_overlap: bool, rna_manifest: Optional[str]) -> pd.DataFrame:
    nrna_dir = normal_root / "normal" / "rna"
    if not nrna_dir.exists():
        return base_df
    rna = _rna_from_manifest(nrna_dir, case_id, rna_manifest)
    if rna is None or rna.empty:
        return base_df
    left_key, right_key = _left_right_keys(ensg_join_mode)
    if synth_overlap:
        need = set(base_df[left_key].dropna().unique()) - set(rna[right_key].dropna().unique())
        if need:
            add = pd.DataFrame({
                "Ensembl_full": [x if ensg_join_mode == "exact" else
                                 (f"{x}.1" if not re.search(r'\.\d+$', str(x)) else str(x)) for x in need],
                "Ensembl_core": list(need), "Count": 100
            })
            rna = pd.concat([rna, add], ignore_index=True)
    merged = base_df.merge(rna[[right_key, "Count"]].rename(columns={"Count": "Normal_RNA_Count"}),
                           how="left", left_on=left_key, right_on=right_key)
    if right_key in merged.columns:
        merged = merged.drop(columns=[right_key])
    logging.info(f"[RNA:NORMAL] {case_id}: matched Normal_RNA_Count for {int(merged['Normal_RNA_Count'].notna().sum())} rows")
    return merged

def integrate_cnv_normal(base_df: pd.DataFrame, normal_root: Path, case_id: str,
                         ensg_join_mode: str, synth_overlap: bool, cn_manifest: Optional[str],
                         ref_dir: Optional[Path]=None) -> pd.DataFrame:
    ncn_dir = normal_root / "normal" / "cn"
    if not ncn_dir.exists():
        ncn_dir = normal_root / "normal" / "cnv"
    if not ncn_dir.exists():
        return base_df
    cn = _cn_from_manifest(ncn_dir, case_id, cn_manifest, ref_dir=ref_dir)
    if cn is None or cn.empty:
        return base_df
    left_key, right_key = _left_right_keys(ensg_join_mode)
    if synth_overlap:
        need = set(base_df[left_key].dropna().unique()) - set(cn[right_key].dropna().unique())
        if need:
            add = pd.DataFrame({
                "Ensembl_full": [x if ensg_join_mode == "exact" else
                                 (f"{x}.1" if not re.search(r'\.\d+$', str(x)) else str(x)) for x in need],
                "Ensembl_core": list(need), "copy_number": 2
            })
            cn = pd.concat([cn, add], ignore_index=True)
    merged = base_df.merge(cn[[right_key, "copy_number"]].rename(columns={"copy_number": "Normal_CNV_Count"}),
                           how="left", left_on=left_key, right_on=right_key)
    if right_key in merged.columns:
        merged = merged.drop(columns=[right_key])
    logging.info(f"[CNV:NORMAL] {case_id}: matched Normal_CNV_Count for {int(merged['Normal_CNV_Count'].notna().sum())} rows")
    return merged

def integrate_ch3_normal(base_df: pd.DataFrame, normal_root: Path, case_id: str, ensg_join_mode: str,
                         ch3_map: Optional[str], ref_dir: Optional[Path],
                         ch3_probe_col: Optional[str], ch3_ensg_col: Optional[str], ch3_symbol_col: Optional[str],
                         ch3_agg: str, ch3_manifest: Optional[str]) -> pd.DataFrame:
    nch3_dir = normal_root / "normal" / "ch3"
    if not nch3_dir.exists():
        return base_df
    ch3 = _ch3_from_manifest(nch3_dir, case_id, ch3_manifest)
    if ch3 is None or ch3.empty:
        return base_df

    m = _load_ch3_map(ch3_map, ref_dir, ch3_probe_col, ch3_ensg_col, ch3_symbol_col)
    if m is None or m.empty:
        logging.warning(f"[CH3:NORMAL] {case_id}: No usable probe→ENSG map; skipping CH3.")
        return base_df

    ch3m = ch3.merge(m, how="left", on="probe").dropna(subset=["ENSGene_core"])
    if ch3m.empty:
        logging.warning(f"[CH3:NORMAL] {case_id}: probe map didn’t cover this case’s probes; skipping.")
        return base_df

    if ch3_agg == "median":
        agg = ch3m.groupby("ENSGene_core", as_index=False)["beta"].median().rename(columns={"beta":"Normal_CH3_Beta"})
    elif ch3_agg == "all":
        stats = ch3m.groupby("ENSGene_core")["beta"].agg(['mean','median','min','max','std','count']).reset_index()
        stats = stats.rename(columns={
            'mean':'Normal_CH3_Beta_mean','median':'Normal_CH3_Beta_median','min':'Normal_CH3_Beta_min',
            'max':'Normal_CH3_Beta_max','std':'Normal_CH3_Beta_std','count':'Normal_CH3_Beta_n','ENSGene_core':'ENSGene_core'
        })
        agg = stats
    else:
        agg = ch3m.groupby("ENSGene_core", as_index=False)["beta"].mean().rename(columns={"beta":"Normal_CH3_Beta"})

    left_key = "ENSGene_core" if ensg_join_mode == "core" else "ENSGene"
    merged = base_df.merge(agg, how="left", left_on=left_key, right_on="ENSGene_core")
    if left_key != "ENSGene_core" and "ENSGene_core" in merged.columns:
        merged = merged.drop(columns=["ENSGene_core"])
    logging.info(f"[CH3:NORMAL] {case_id}: matched Normal_CH3_Beta for {int(merged.filter(like='Normal_CH3_Beta').notna().sum().sum())} rows")
    return merged

def integrate_protein_normal(base_df: pd.DataFrame, normal_root: Path, case_id: str,
                             protein_policy: str = "best") -> pd.DataFrame:
    nprot_dir = normal_root / "normal" / "protein"
    p = _try_case_file_in_dir(nprot_dir, case_id)
    if p is None:
        return base_df
    prot = _read_any_table(p)  # auto sep by extension
    prot = prot.drop_duplicates()
    if "NP" not in prot.columns:
        logging.warning(f"[PROTEIN:NORMAL] {case_id}: protein file lacks NP; skipping.")
        return base_df
    if protein_policy == "best":
        prot = _select_best_protein_rows(prot)
    elif protein_policy == "first":
        if "SEQ" in prot.columns:
            prot = prot[prot["SEQ"].astype(str).str.len() > 3]
        prot = prot.drop_duplicates(subset=["NP"], keep="first")
    else:
        if "SEQ" in prot.columns:
            prot = prot[prot["SEQ"].astype(str).str.len() > 3]

    prot = prot.rename(columns={
        "SEQ": "Normal_SEQ",
        "EV":  "Normal_EV",
        "INT": "Normal_INT",
    })
    df = base_df.copy()
    # Handle semicolon-joined NP lists by exploding
    if "NP" in df.columns and df["NP"].astype(str).str.contains(";").any():
        df = df.copy()
        df["_orig_idx"] = range(len(df))
        df["NP"] = df["NP"].astype(str).apply(lambda x: [t.strip() for t in x.split(";") if t.strip()] if isinstance(x, str) else [])
        df = df.explode("NP", ignore_index=True)
        merged = df.merge(prot[["NP","Normal_SEQ","Normal_EV","Normal_INT"]], how="left", on="NP")
        # keep best per original row (lowest EV, longest seq)
        if "Normal_EV" in merged.columns:
            merged["_ev"] = pd.to_numeric(merged["Normal_EV"], errors="coerce")
        else:
            merged["_ev"] = float("inf")
        if "Normal_SEQ" in merged.columns:
            merged["_len"] = merged["Normal_SEQ"].astype(str).str.len()
        else:
            merged["_len"] = -1
        merged = merged.sort_values(["_orig_idx", "_ev", "_len"], ascending=[True, True, False])
        merged = merged.drop_duplicates(subset=["_orig_idx"], keep="first")
        merged = merged.drop(columns=["_ev","_len","_orig_idx"])
        df = merged
    else:
        df["NP"] = df.get("NP", pd.Series(index=df.index, dtype=object)).astype(str).str.strip()
        df = df.merge(prot[["NP","Normal_SEQ","Normal_EV","Normal_INT"]], how="left", on="NP")
    logging.info(f"[PROTEIN:NORMAL] {case_id}: matched Normal_SEQ for {int(df['Normal_SEQ'].notna().sum()) if 'Normal_SEQ' in df.columns else 0} rows")
    return df

# ------------------------- Per-case worker --------------------------

def process_one_case(case_id: str,
                     dna_dir: str, rna_dir: str, ch3_dir: str, cn_dir: str, protein_dir: str,
                     out_dir: str, ref_dir: Optional[str],
                     step: str, dna_rows: str, dedup_level: str, dedup_key_cols: Optional[List[str]],
                     emit_qc: bool, ensg_join_mode: str, synth_overlap: bool,
                     ch3_map: Optional[str], ch3_probe_col: Optional[str], ch3_ensg_col: Optional[str],
                     ch3_symbol_col: Optional[str], ch3_agg: str,
                     skip_rna: bool, skip_ch3: bool, skip_protein: bool, skip_cn: bool,
                     rna_manifest: Optional[str], ch3_manifest: Optional[str], cn_manifest: Optional[str],
                     keep_join_cols: bool,
                     normal_dir: Optional[str] = None,
                     protein_policy: str = "best"
                     ) -> Tuple[str, Optional[str], dict]:
    try:
        out_dir_p = Path(out_dir)
        ref_p = Path(ref_dir) if ref_dir else None
        dna_dir_p = Path(dna_dir); rna_dir_p = Path(rna_dir)
        ch3_dir_p = Path(ch3_dir); cn_dir_p = Path(cn_dir)
        protein_dir_p = Path(protein_dir)

        # Build DNA if requested or missing
        dna_out_dir = out_dir_p / "dna"; _ensure_dir(dna_out_dir)
        base_fp = dna_out_dir / f"{case_id}.csv"
        if step in ("dna", "all") or not base_fp.exists():
            load_and_enhance_dna(case_id, dna_dir_p, out_dir_p, ref_p, dna_rows)
        base = _read_any_table(base_fp)

        # Decide which integrations to run
        chain = (step == "dna")
        do_rna = ((step in ("rna", "all")) or chain) and (not skip_rna)
        do_ch3 = ((step in ("ch3", "all")) or chain) and (not skip_ch3)
        do_pro = ((step in ("protein", "all")) or chain) and (not skip_protein)
        do_cnv = ((step in ("cnv", "all")) or chain) and (not skip_cn)

        # ---- Tumor merges
        if do_rna:
            base = integrate_rna(base, rna_dir_p, case_id, ensg_join_mode, synth_overlap, rna_manifest)
        if do_ch3:
            base = integrate_ch3(base, ch3_dir_p, case_id, ensg_join_mode, ch3_map, ref_p,
                                 ch3_probe_col, ch3_ensg_col, ch3_symbol_col, ch3_agg, ch3_manifest)
        if do_pro:
            base = integrate_protein(base, protein_dir_p, case_id, synth_overlap, protein_policy=protein_policy)
        if do_cnv:
            base = integrate_cnv(base, cn_dir_p, case_id, ensg_join_mode, synth_overlap, cn_manifest, ref_dir=ref_p)

        # Standardize names
        rename_map = {"Count": "RNA_Count", "beta_val": "CH3_Beta", "copy_number": "CNV_Count"}
        base = base.rename(columns={k: v for k, v in rename_map.items() if k in base.columns})

        # Dedup before normal
        base = _final_dedup(base, dedup_level, dedup_key_cols)

        # ---- Normal merges (actual integration)
        if normal_dir:
            normal_root = Path(normal_dir)
            if do_rna:
                base = integrate_rna_normal(base, normal_root, case_id, ensg_join_mode, synth_overlap, rna_manifest)
            if do_ch3:
                base = integrate_ch3_normal(base, normal_root, case_id, ensg_join_mode,
                                            ch3_map, ref_p, ch3_probe_col, ch3_ensg_col, ch3_symbol_col,
                                            ch3_agg, ch3_manifest)
            if do_pro:
                base = integrate_protein_normal(base, normal_root, case_id, protein_policy=protein_policy)
            if do_cnv:
                base = integrate_cnv_normal(base, normal_root, case_id, ensg_join_mode, synth_overlap, cn_manifest, ref_dir=ref_p)

        
        # ---- Metrics for summary (compute before dropping join columns)
        metrics = {
            "case_id": case_id,
            "dna_total_rows": int(len(base)),
            "dna_ensg_rows": int(base["ENSGene"].notna().sum()) if "ENSGene" in base.columns else (
                int(base["ENSGene_core"].notna().sum()) if "ENSGene_core" in base.columns else None
            ),
            "tumor_rna_mapped": int(base["RNA_Count"].notna().sum()) if "RNA_Count" in base.columns else 0,
            "tumor_cn_mapped": int(base["CNV_Count"].notna().sum()) if "CNV_Count" in base.columns else (
                int(base["copy_number"].notna().sum()) if "copy_number" in base.columns else 0
            ),
            "tumor_protein_mapped": int(base["SEQ"].notna().sum()) if "SEQ" in base.columns else (
                int(base["Protein_SEQ"].notna().sum()) if "Protein_SEQ" in base.columns else 0
            ),
            "tumor_ch3_mapped": int(base["CH3_Beta"].notna().sum()) if "CH3_Beta" in base.columns else 0,
            "normal_rna_mapped": int(base["Normal_RNA_Count"].notna().sum()) if "Normal_RNA_Count" in base.columns else 0,
            "normal_cn_mapped": int(base["Normal_CNV_Count"].notna().sum()) if "Normal_CNV_Count" in base.columns else (
                int(base["Normal_copy_number"].notna().sum()) if "Normal_copy_number" in base.columns else 0
            ),
            "normal_protein_mapped": int(base["Normal_SEQ"].notna().sum()) if "Normal_SEQ" in base.columns else 0,
            "normal_ch3_mapped": int(base["Normal_CH3_Beta"].notna().sum()) if "Normal_CH3_Beta" in base.columns else 0,
            "errors": ""
        }

        # ---- Soft cleanup (drop join/reference columns) unless asked to keep them
        if not keep_join_cols:
            linker_exact = {
                "ENSGene", "ENSGene_core", "Ensembl_core", "Ensembl_full",
                "symbol", "probe", "To_core", "From", "Gene_Name", "UniProt", "NP", "CpG_Probes"
            }
            linker_suffixes = tuple(["_core", "_probe"])
            to_drop = [c for c in base.columns
                       if c in linker_exact or any(c.endswith(sfx) for sfx in linker_suffixes)]
            to_drop = [c for c in to_drop if c not in {"Gene"}]  # keep readable symbol if present
            if to_drop:
                base = base.drop(columns=sorted(set(to_drop)), errors="ignore")

            # Reorder NP and CNV_Count to appear just before SEQ if present (tumor)
            if "SEQ" in base.columns:
                cols = list(base.columns)
                if "NP" in cols:
                    cols.remove("NP"); seq_idx = cols.index("SEQ"); cols.insert(seq_idx, "NP")
                if "CNV_Count" in cols:
                    cols.remove("CNV_Count"); seq_idx = cols.index("SEQ"); cols.insert(seq_idx, "CNV_Count")
                base = base[cols]

        # Write
        out_fp = out_dir_p / f"{case_id}_integrated.csv"
        base.to_csv(out_fp, index=False)

        if emit_qc:
            qc = {
                "case_id": case_id,
                "n_rows_out": len(base),
                "rna_nonnull": int(base.filter(regex=r"^RNA_Count$").notna().sum().sum()),
                "cnv_nonnull": int(base.filter(regex=r"^CNV_Count$").notna().sum().sum()),
                "ch3_nonnull": int(base.filter(regex=r"^CH3_Beta$|^CH3_Beta_").notna().sum().sum()),
                "prot_nonnull": int(base.filter(regex=r"^SEQ$").notna().sum().sum()),
                "normal_rna_nonnull": int(base.filter(regex=r"^Normal_RNA_Count$").notna().sum().sum()),
                "normal_cnv_nonnull": int(base.filter(regex=r"^Normal_CNV_Count$").notna().sum().sum()),
                "normal_ch3_nonnull": int(base.filter(regex=r"^Normal_CH3_Beta$|^Normal_CH3_Beta_").notna().sum().sum()),
                "normal_prot_nonnull": int(base.filter(regex=r"^Normal_SEQ$").notna().sum().sum()),
            }
            # per-case QC disabled (no qc dir created)

        _emit_qc_counts(case_id, base)

        return (case_id, None, metrics)
    except Exception as e:
        return (case_id, str(e), {"case_id": case_id, "errors": str(e)})

_DNA_RE = re.compile(r'\[DNA\]\s+(\S+):\s+ENSG extracted for\s+(\d+)\s*/\s*(\d+)\s+rows')
_DNA_FILTER_RE = re.compile(r'\[DNA\]\s+(\S+):\s+filtered to ENSG-only rows:\s+(\d+)\s*/\s*(\d+)\s*')

_RNA_T_RE = re.compile(r'\[RNA\]\s+(\S+):\s+matched\s+(?:Count|RNA_Count)\s+for\s+(\d+)\s+rows')
_CNV_T_RE = re.compile(r'\[CNV\]\s+(\S+):\s+matched\s+(?:copy_number|CNV_Count)\s+(?:for|in)\s+(\d+)\s+rows')
_PRO_T_RE = re.compile(r'\[PROTEIN\]\s+(\S+):\s+matched\s+NP\s+for\s+(\d+)\s+rows')
_CH3_T_RE = re.compile(r'\[CH3\]\s+(\S+):\s+matched\s+(?:beta_val|CH3_Beta)\s+for\s+(\d+)\s+rows')

_RNA_N_RE = re.compile(r'\[RNA:NORMAL\]\s+(\S+):\s+matched\s+Normal_RNA_Count\s+for\s+(\d+)\s+rows')
_CNV_N_RE = re.compile(r'\[CNV:NORMAL\]\s+(\S+):\s+matched\s+Normal_CNV_Count\s+for\s+(\d+)\s+rows')
_PRO_N_RE = re.compile(r'\[PROTEIN:NORMAL\]\s+(\S+):\s+matched\s+Normal_SEQ\s+for\s+(\d+)\s+rows')
_CH3_N_RE = re.compile(r'\[CH3:NORMAL\]\s+(\S+):\s+matched\s+Normal_CH3_Beta\s+for\s+(\d+)\s+rows')

_ERR_RE    = re.compile(r'ERROR \[(\S+)\]:\s+(.*)')

def _summarize_log(log_path: Path) -> pd.DataFrame:
    cases = {}

    def ensure(k):
        if k not in cases:
            cases[k] = dict(
                case_id=k, dna_total_rows=None, dna_ensg_rows=None,
                tumor_rna_mapped=0, tumor_cn_mapped=0, tumor_protein_mapped=0, tumor_ch3_mapped=0,
                normal_rna_mapped=0, normal_cn_mapped=0, normal_protein_mapped=0, normal_ch3_mapped=0,
                errors=[]
            )
        return cases[k]

    with log_path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            for regex, setter in (
                (_DNA_RE,      lambda m: (ensure(m[1]).update(dna_ensg_rows=int(m[2]), dna_total_rows=int(m[3])))),
                (_DNA_FILTER_RE, lambda m: (ensure(m[1]).update(dna_ensg_rows=ensure(m[1]).get("dna_ensg_rows") or int(m[2]),
                                                     dna_total_rows=ensure(m[1]).get("dna_total_rows") or int(m[3])))),
                (_RNA_T_RE,    lambda m: (ensure(m[1]).update(tumor_rna_mapped=int(m[2])))),
                (_CNV_T_RE,    lambda m: (ensure(m[1]).update(tumor_cn_mapped=int(m[2])))),
                (_PRO_T_RE,    lambda m: (ensure(m[1]).update(tumor_protein_mapped=int(m[2])))),
                (_CH3_T_RE,    lambda m: (ensure(m[1]).update(tumor_ch3_mapped=int(m[2])))),
                (_RNA_N_RE,    lambda m: (ensure(m[1]).update(normal_rna_mapped=int(m[2])))),
                (_CNV_N_RE,    lambda m: (ensure(m[1]).update(normal_cn_mapped=int(m[2])))),
                (_PRO_N_RE,    lambda m: (ensure(m[1]).update(normal_protein_mapped=int(m[2])))),
                (_CH3_N_RE,    lambda m: (ensure(m[1]).update(normal_ch3_mapped=int(m[2])))),
            ):
                m = regex.search(line)
                if m:
                    setter(m)
                    break
            else:
                m = _ERR_RE.search(line)
                if m:
                    ensure(m[1])["errors"].append(m[2].strip())

    df = pd.DataFrame(cases.values())
    if df.empty:
        df = pd.DataFrame(columns=[
            'case_id','dna_total_rows','dna_ensg_rows',
            'tumor_rna_mapped','tumor_cn_mapped','tumor_protein_mapped','tumor_ch3_mapped',
            'normal_rna_mapped','normal_cn_mapped','normal_protein_mapped','normal_ch3_mapped','errors']
        )
        return df
    df = df.sort_values('case_id').reset_index(drop=True)
    if "errors" in df.columns:
        df["errors"] = df["errors"].map(lambda x: "; ".join(x) if isinstance(x, list) and x else "")
    return df

def write_run_summary(log_path: Path, out_path: Path, fmt: str = None):
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df = _summarize_log(log_path)
    if fmt is None:
        fmt = "csv" if out_path.suffix.lower() == ".csv" else "tsv"
    if fmt == "csv":
        df.to_csv(out_path, index=False)
    else:
        df.to_csv(out_path, sep="\t", index=False)
    logging.info(f"[QC] Wrote run summary: {out_path} ({len(df)} cases)")

def write_qc_cohort_from_files(qc_dir: Path, out_path: Path):
    """Combine per-case QC TSVs into a single CSV/TSV."""
    fps = sorted(qc_dir.glob("*_qc.tsv"))
    if not fps:
        logging.warning(f"[QC] No per-case QC files found under {qc_dir} to combine.")
        return
    dfs = []
    for fp in fps:
        try:
            df = pd.read_csv(fp, sep="\t", dtype=str)
            dfs.append(df)
        except Exception as e:
            logging.error(f"[QC] Failed to read {fp}: {e}")
    if not dfs:
        logging.warning(f"[QC] No readable QC files under {qc_dir}")
        return
    big = pd.concat(dfs, ignore_index=True)
    # Cast numeric-looking columns
    for col in big.columns:
        if col != "case_id":
            big[col] = pd.to_numeric(big[col], errors="coerce")
    # reorder columns
    order = ['case_id','n_rows_out','rna_nonnull','cnv_nonnull','prot_nonnull','ch3_nonnull',
             'normal_rna_nonnull','normal_cnv_nonnull','normal_prot_nonnull','normal_ch3_nonnull']
    cols = [c for c in order if c in big.columns] + [c for c in big.columns if c not in order]
    big = big[cols]
    out_path.parent.mkdir(parents=True, exist_ok=True)
    if out_path.suffix.lower() == ".csv":
        big.to_csv(out_path, index=False)
    else:
        big.to_csv(out_path, sep="\t", index=False)
    logging.info(f"[QC] Wrote cohort QC: {out_path} ({len(big)} cases)")

# ------------------------------- CLI --------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Multiomic integration: DNA core + optional RNA/CH3/CNV/Protein (+normal).")
    p.add_argument("--in_dir", required=True, help="Project root directory")
    p.add_argument("--manifest", required=True, help="Case list: first column contains case IDs")
    p.add_argument("--add-normal-dir", dest="add_normal_dir", default=None,
                   help="Root containing normal/{rna,ch3,cn,protein} to append Normal_* columns.")
    p.add_argument("--normal-dir", dest="add_normal_dir", default=None,
                   help="Alias of --add-normal-dir.")  # alias
    p.add_argument("--out_dir", required=True, help="Output directory")
    p.add_argument("--ref_dir", default=None, help="Reference directory containing mapping tables")
    p.add_argument("--ref_zip", default=None, help="Zip archive with reference tables (extracts to temp)")
    p.add_argument("--step", default="all", choices=["all", "dna", "rna", "ch3", "cnv", "protein"], help="Pipeline step(s)")
    p.add_argument("--jobs", type=int, default=os.cpu_count() or 2, help="Parallel workers")
    p.add_argument("--case-id", default=None, help="Process a single case-id instead of all in --manifest")

    # DNA handling / output controls
    p.add_argument("--dna_rows", default="ensg_only", choices=["ensg_only", "all"],
                   help="Rows to keep from DNA before joins")
    p.add_argument("--dedup_level", default="row", choices=["none", "row", "key"],
                   help="Final output de-duplication level")
    p.add_argument("--dedup_key", default="CHROM,POS,REF,ALT,ENSGene,ENSGene_core",
                   help="Columns to use when --dedup_level=key")
    p.add_argument("--emit_qc", action="store_true", help="Write per-case QC to out_dir/qc/*.tsv")

    # Join behavior
    p.add_argument("--ensg_join_mode", default="core", choices=["core", "exact"],
                   help="Join by versionless ENSG (core) or exact IDs")
    p.add_argument("--synthesize_overlap_for_tests", action="store_true",
                   help="Fill missing RNA/CNV/Protein rows for NP/ENSG present in DNA to ease tests")

    # CH3 options
    p.add_argument("--ch3_map", default=None, help="Probe→gene map (CSV/TSV). If not set, try reference.zip.")
    p.add_argument("--ch3_probe_col", default=None, help="Column name for probe IDs in --ch3_map")
    p.add_argument("--ch3_ensg_col", default=None, help="Column name for ENSG IDs in --ch3_map")
    p.add_argument("--ch3_symbol_col", default=None, help="Column name for gene symbols in --ch3_map")
    p.add_argument("--ch3_agg", default="mean", choices=["mean", "median", "all"], help="Aggregate beta per gene")

    # Protein policy
    p.add_argument("--protein_policy", default="best", choices=["best", "first", "all"],
                   help="How to select protein rows per NP: best=lowest EV then longest SEQ; first=first per NP; all=keep all")

    # Skip flags
    p.add_argument("--skip-rna", action="store_true", help="Skip RNA integration when chaining from --step dna or all.")
    p.add_argument("--skip-ch3", action="store_true", help="Skip CH3 integration when chaining from --step dna or all.")
    p.add_argument("--skip-protein", action="store_true", help="Skip protein integration when chaining from --step dna or all.")
    p.add_argument("--skip-cn", action="store_true", help="Skip CNV integration when chaining from --step dna or all.")

    # Manifests for raw data (optional; will auto-detect if omitted)
    p.add_argument("--rna-manifest", dest="rna_manifest", default=None,
                   help="Path to RNA manifest (case→file). Defaults to autodetect in <input_rna_dir>/")
    p.add_argument("--ch3-manifest", dest="ch3_manifest", default=None,
                   help="Path to CH3 manifest (case→file). Defaults to autodetect in <input_ch3_dir>/")
    p.add_argument("--cn-manifest", dest="cn_manifest", default=None,
                   help="Path to CNV manifest (case→file). Defaults to autodetect in <input_cn_dir>/")

    # Inputs for raw data
    p.add_argument("--input_dna_dir", default=None,
                   help="Directory containing DNA per-case files (default: <in_dir>/dna)")
    p.add_argument("--input_rna_dir", default=None,
                   help="Directory containing RNA files/in_dirs (default: <in_dir>/rna)")
    p.add_argument("--input_ch3_dir", default=None,
                   help="Directory containing methylation files/in_dirs (default: <in_dir>/ch3)")
    p.add_argument("--input_cn_dir", default=None,
                   help="Directory containing CNV files/in_dirs (default: <in_dir>/cnv or <in_dir>/cn)")
    p.add_argument("--input_protein_dir", default=None,
                   help="Directory containing protein per-case files (default: <in_dir>/protein)")
    p.add_argument("--log-file", default=None,
                   help="Write run logs to this file (default: <out_dir>/metadata/multiomic_run.log)")
    p.add_argument("--summary-out", default=None,
                   help="Path for the run summary table (default: <out_dir>/metadata/qc_summary.tsv)")
    p.add_argument("--summary-format", choices=["tsv", "csv"], default=None,
                   help="Summary format (default inferred from --summary-out; else TSV)")

    # Optional: combined per-case QC into cohort CSV/TSV
    p.add_argument("--qc-cohort-out", default=None,
                   help="Write a single cohort QC file by combining per-case qc/*.tsv (default: <out_dir>/qc_cohort.csv)")

    # Soft cleanup control
    p.add_argument("--keep-join-cols", action="store_true",
                   help="Keep intermediate join keys (e.g., ENSGene_core, Ensembl_core, probe). "
                        "By default they are dropped and integrated columns are renamed.")
    return p.parse_args()

def main() -> None:
    args = parse_args()
    out_dir = Path(args.out_dir); _ensure_dir(out_dir)

    log_path = Path(args.log_file) if args.log_file else Path(args.out_dir) / "metadata" / "multiomic_run.log"
    log_path.parent.mkdir(parents=True, exist_ok=True)

    fmt = "[multiomic_integration] %(levelname)s %(message)s"
    handlers = [logging.StreamHandler(sys.stdout), logging.FileHandler(log_path, mode="w", encoding="utf-8")]
    logging.basicConfig(level=logging.INFO, format=fmt, handlers=handlers)

    logging.info(f"[AUTO] Logging to {log_path}")

    # Resolve reference directory
    ref_dir_p = resolve_reference_dir(args.ref_dir, args.ref_zip)
    if ref_dir_p is not None:
        logging.info(f"Using reference dir for maps: {ref_dir_p}")

    # Case selection
    case_ids = [args.case_id] if args.case_id else load_case_ids(args.manifest)

    # Dedup key columns
    dedup_key_cols = [c.strip() for c in args.dedup_key.split(",") if c.strip()] if args.dedup_level == "key" else None

    in_dir_p = Path(args.in_dir)

    # Input dirs (allow overrides)
    dna_dir = Path(args.input_dna_dir) if args.input_dna_dir else (in_dir_p / "dna")
    rna_dir = Path(args.input_rna_dir) if args.input_rna_dir else (in_dir_p / "rna")
    ch3_dir = Path(args.input_ch3_dir) if args.input_ch3_dir else (in_dir_p / "ch3")
    protein_dir = Path(args.input_protein_dir) if args.input_protein_dir else (in_dir_p / "protein")

    if args.input_cn_dir:
        cn_dir = Path(args.input_cn_dir)
    else:
        cn_dir = in_dir_p / "cnv"
        if not cn_dir.exists():
            cn_dir = in_dir_p / "cn"

    logging.info(f"Dirs: dna={dna_dir} rna={rna_dir} ch3={ch3_dir} cn={cn_dir} protein={protein_dir}")
    logging.info(f"N = {len(case_ids)} samples; step(s)={args.step}")

    with ProcessPoolExecutor(max_workers=args.jobs) as ex:
        futs = [ex.submit(
            process_one_case, cid, str(dna_dir), str(rna_dir), str(ch3_dir), str(cn_dir), str(protein_dir),
            str(out_dir), str(ref_dir_p) if ref_dir_p else None,
            args.step, args.dna_rows, args.dedup_level, dedup_key_cols, args.emit_qc,
            args.ensg_join_mode, args.synthesize_overlap_for_tests,
            args.ch3_map, args.ch3_probe_col, args.ch3_ensg_col, args.ch3_symbol_col, args.ch3_agg,
            args.skip_rna, args.skip_ch3, args.skip_protein, args.skip_cn,
            args.rna_manifest, args.ch3_manifest, args.cn_manifest,
            args.keep_join_cols,
            args.add_normal_dir,
            args.protein_policy
        ) for cid in case_ids]

        for cid in case_ids:
            logging.info(f"Processing {cid} [step={args.step}]")
        
        metrics_rows = []
        for f in as_completed(futs):
            cid, err, met = f.result()
            if err:
                logging.error(f"[{cid}]: {err}")
            if met:
                metrics_rows.append(met)
        logging.info("Done.")

        # Write summary directly from metrics_rows
        import pandas as pd
        summary_path = Path(args.summary_out) if args.summary_out else Path(args.out_dir) / "metadata" / "qc_summary.tsv"
        fmt = args.summary_format
        summary_path.parent.mkdir(parents=True, exist_ok=True)
        if metrics_rows:
            df_sum = pd.DataFrame(metrics_rows).sort_values("case_id")
        else:
            df_sum = pd.DataFrame(columns=[
                "case_id","dna_total_rows","dna_ensg_rows",
                "tumor_rna_mapped","tumor_cn_mapped","tumor_protein_mapped","tumor_ch3_mapped",
                "normal_rna_mapped","normal_cn_mapped","normal_protein_mapped","normal_ch3_mapped","errors"
            ])
        if fmt == "csv" or (fmt is None and summary_path.suffix.lower()==".csv"):
            df_sum.to_csv(summary_path, index=False)
        else:
            df_sum.to_csv(summary_path, sep="\t", index=False)
        logging.info(f"[QC] Wrote run summary: {summary_path} ({len(df_sum)} cases)")

        # Cohort QC from summary (guard existence)
        qc_out = Path(args.qc_cohort_out) if args.qc_cohort_out else (out_dir / "metadata" / "qc_cohort.csv")
        if summary_path.exists():
            try:
                if summary_path.suffix.lower() == ".csv":
                    _sum = pd.read_csv(summary_path, dtype=str)
                else:
                    _sum = pd.read_csv(summary_path, sep="\t", dtype=str)
                qc_out.parent.mkdir(parents=True, exist_ok=True)
                _sum.to_csv(qc_out, index=False)
                logging.info(f"[QC] Wrote cohort QC from summary: {qc_out} ({len(_sum)} cases)")
            except Exception as e:
                logging.error(f"[QC] Failed to write cohort QC from summary: {e}")
        else:
            logging.error(f"[QC] Summary file not found: {summary_path}")

        # --- Cleanup: remove empty out_dir/qc if it exists and is empty ---
        try:
            qcd = out_dir / "qc"
            if qcd.exists():
                # if directory is empty, remove it
                if not any(qcd.iterdir()):
                    qcd.rmdir()
                    logging.info(f"[QC] Removed empty directory: {qcd}")
        except Exception as _e:
            logging.warning(f"[QC] Could not remove qc dir: { _e }")
        return

if __name__ == "__main__":
    main()
