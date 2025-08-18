#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MAKE MANIFEST

This script organizes metafiles (manifests) that will be used to integrate multiomic information for calling during multiomic integration. Outputs are tab-separated TSVs.

Dependencies: pandas

Project layout expected (for --build main)

PROJECT_ROOT/
   dna/        (optional)         case-id.csv                       #VEP-annotated mutect files
               gdc-dna.tsv        + dna/{File ID}/{File Name}
   rna/        gdc-rna.tsv        + rna/{File ID}/{File Name}       #RNA-seq count data
   ch3/        gdc-ch3.tsv        + ch3/{File ID}/{File Name}       #Methylation data in GDC format
   cn/         gdc-cn.tsv         + cn/{File ID}/{File Name}        #Copy Number data in GDC format
   protein/    (optional)         case-id.csv                       #PSM reads (NP, SEQ, EV, INT)

Usage 
(Recommended):
   python make_manifest.py \
        --build-main \
        --project_root <project directory> \
        --out_dir <output directory>

(Using existing unified manifest):
   python make_manifest.py \
        --emit_ref \
        --main-manifest <main-manifest.tsv>
        --project_root <project directory> \
        --out_dir <output directory>

(To make single manifest from gdc-manifest):
   python make_manifest.py \
        --gdc-manifest <manifest file> \
        --out_dir <output directory> \
        --type <rna|ch3|cn> \
        --folder-id-col <# (e.g., 0)> \
        --file-name-col <# (e.g., 1)> \
        --case-id-col <# (e.g, 5)> 

Arguments:
(Recommended Use)
   --build-main        Flag that creates unified manifest with paths constructed like {PROJECT_ROOT}/{modality}/{File ID}/{File Name}
   --project_root      Directory that has all project-related files (dna, rna, ch3, cn, protein)
   --main_out          Directory that will contain all output files

(General)
   --log-level L       Default: INFO, can be set to DEBUG for more verbose logging
   
(Single Use)   
   --gdc-manifest      Flag to use if opting for creation of single manifest as opposed to build-main
   --out_dir           Directory that will contain output file
   --type              Used for naming output manifest (e.g., rna_manifest.tsv)
   --folder-id-column  Column number that identifies path (downloaded code path; e.g, 123abc)
   --file-name-col     Column number that identifies file name (downloaded code + file name; e.g., 987zyx.wxs.MuTect2.somatic_annotation.vcf.gz)
   --case-id-col       Column number that identifies case-id (e.g., C3N-001)
"""

from __future__ import annotations
import argparse
import logging
from pathlib import Path
import re
import pandas as pd

logger = logging.getLogger("make_manifest")

# ---------- Helpers ----------

def _norm_col(s: str) -> str:
    """Normalize column keys: strip, lowercase, collapse whitespace/underscores to hyphen."""
    return re.sub(r"[\s_]+", "-", str(s).strip().lower())

def _read_any(path: Path) -> pd.DataFrame:
    """Read CSV/TSV with autodetected delimiter and normalize headers."""
    df = pd.read_csv(path, sep=None, engine="python")
    df.columns = [_norm_col(c) for c in df.columns]
    return df

def _read_gdc(path: Path) -> pd.DataFrame:
    """Read a GDC TSV with canonical headers (no normalization)."""
    return pd.read_csv(path, sep="\t")

def _write_tsv(df: pd.DataFrame, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    rows, cols = getattr(df, "shape", (0, 0))
    logger.info("Writing TSV to %s (rows=%d, cols=%d)", out_path, rows, cols)
    df.to_csv(out_path, sep="\t", index=False)
    logger.info("Wrote: %s", out_path)

def _norm_case_id(val: str) -> str:
    """Return the first token before a comma; strip quotes/spaces."""
    s = str(val).strip().strip('"').strip("'")
    if "," in s:
        s = s.split(",")[0].strip()
    return s

def _infer_case_column(cols) -> str | None:
    """Infer a Case ID column name using common variants."""
    norm = {c.lower().replace("_", " ").replace("-", " ").strip(): c for c in cols}
    for k in ("case id", "caseid"):
        if k in norm:
            return norm[k]
    for c in cols:
        lc = c.lower()
        if lc in ("case-id", "case_id", "caseid", "case id"):
            return c
    return None

# -------------------- RNA --------------------

def _infer_rna_path_column(cols) -> str | None:
    preferred = ("rna_path", "rna_file", "rna", "rna_csv")
    for p in preferred:
        for c in cols:
            if c.lower() == p:
                return c
    for c in cols:
        lc = c.lower()
        if "rna" in lc and ("path" in lc or "file" in lc or "csv" in lc):
            return c
    return None

def _write_rna_manifest_from_main(main_manifest_path: Path, rna_manifest_path: Path) -> int:
    df = pd.read_csv(main_manifest_path, sep="\t")
    case_col = _infer_case_column(df.columns)
    if not case_col:
        raise ValueError("Could not infer Case ID column for RNA manifest")
    rna_col = _infer_rna_path_column(df.columns)
    if not rna_col:
        logger.info("No RNA path column found; skipping rna_manifest.tsv")
        return 0
    out = df[[case_col, rna_col]].copy()
    out.columns = ["case-id", "path"]
    out = out[out["path"].astype(str).str.strip() != ""]
    if out.empty:
        logger.info("RNA path column empty; skipping rna_manifest.tsv")
        return 0
    _write_tsv(out, rna_manifest_path)
    return len(out)

# -------------------- CH3 --------------------

def _infer_ch3_path_column(cols) -> str | None:
    preferred = ("ch3_path", "ch3_file", "ch3", "ch3_csv")
    for p in preferred:
        for c in cols:
            if c.lower() == p:
                return c
    for c in cols:
        lc = c.lower()
        if "ch3" in lc and ("path" in lc or "file" in lc or "csv" in lc):
            return c
    return None

def _write_ch3_manifest_from_main(main_manifest_path: Path, ch3_manifest_path: Path) -> int:
    df = pd.read_csv(main_manifest_path, sep="\t")
    case_col = _infer_case_column(df.columns)
    if not case_col:
        raise ValueError("Could not infer Case ID column for CH3 manifest")
    ch3_col = _infer_ch3_path_column(df.columns)
    if not ch3_col:
        logger.info("No CH3 path column found; skipping ch3_manifest.tsv")
        return 0
    out = df[[case_col, ch3_col]].copy()
    out.columns = ["case-id", "path"]
    out = out[out["path"].astype(str).str.strip() != ""]
    if out.empty:
        logger.info("CH3 path column empty; skipping ch3_manifest.tsv")
        return 0
    _write_tsv(out, ch3_manifest_path)
    return len(out)

# -------------------- CN --------------------

def _infer_cn_path_column(cols) -> str | None:
    preferred = ("cn_path", "cn_file", "cn", "cn_csv", "copy_number")
    for p in preferred:
        for c in cols:
            if c.lower() == p:
                return c
    for c in cols:
        lc = c.lower()
        if ("cn" in lc or "copy_number" in lc) and ("path" in lc or "file" in lc or "csv" in lc):
            return c
    return None

def _write_cn_manifest_from_main(main_manifest_path: Path, cn_manifest_path: Path) -> int:
    df = pd.read_csv(main_manifest_path, sep="\t")
    case_col = _infer_case_column(df.columns)
    if not case_col:
        raise ValueError("Could not infer Case ID column for CN manifest")
    cn_col = _infer_cn_path_column(df.columns)
    if not cn_col:
        logger.info("No CN path column found; skipping cn_manifest.tsv")
        return 0
    out = df[[case_col, cn_col]].copy()
    out.columns = ["case-id", "path"]
    out = out[out["path"].astype(str).str.strip() != ""]
    if out.empty:
        logger.info("CN path column empty; skipping cn_manifest.tsv")
        return 0
    _write_tsv(out, cn_manifest_path)
    return len(out)

# ---------- Build a unified main-manifest.tsv ----------

def build_main_manifest(project_root: Path, out_path: Path, prefer_sample_type: str = "Primary Tumor") -> Path:
    """
    Build a combined manifest with columns: Case-ID, DNA, RNA, CH3, CN, Protein (if available).

    Looks for per-modality GDC TSVs at:
      rna/gdc-rna.tsv
      ch3/gdc-ch3.tsv
      cn/gdc-cn.tsv  (or cnv/gdc-cn.tsv)
      vep/gdc-vep.tsv  (preferred for DNA) or dna/gdc-dna.tsv

    Paths in the unified file look like:
      {project_root}/{modality}/{File ID}/{File Name}

    Protein files are picked by filename from: protein/{Case-ID}.csv (or any *.csv containing the Case-ID).
    """
    project_root = Path(project_root).resolve()
    logger.info("Building main manifest from root: %s", project_root)
   # inside build_main_manifest(...)
def _first_existing(paths):
    for p in paths:
        if p and Path(p).exists():
            return Path(p)
    return None

def _rglob_first(root: Path, pattern: str):
    hits = list(root.rglob(pattern))
    return hits[0] if hits else None

logger.info("Looking for modality GDC tables under %s", project_root)

candidates_map = {
    "rna": [
        project_root / "rna" / "gdc-rna.tsv",
        _rglob_first(project_root, "rna/gdc-rna.tsv"),
        _rglob_first(project_root, "gdc-rna.tsv"),
    ],
    "ch3": [
        project_root / "ch3" / "gdc-ch3.tsv",
        _rglob_first(project_root, "ch3/gdc-ch3.tsv"),
        _rglob_first(project_root, "gdc-ch3.tsv"),
    ],
    "cn": [
        project_root / "cn" / "gdc-cn.tsv",
        project_root / "cnv" / "gdc-cn.tsv",
        _rglob_first(project_root, "cn/gdc-cn.tsv"),
        _rglob_first(project_root, "cnv/gdc-cn.tsv"),
        _rglob_first(project_root, "gdc-cn.tsv"),
    ],
    "dna": [
        project_root / "vep" / "gdc-vep.tsv",   # preferred
        project_root / "dna" / "gdc-vep.tsv",   # <-- accept dna/gdc-vep.tsv too
        project_root / "dna" / "gdc-dna.tsv",
        _rglob_first(project_root, "vep/gdc-vep.tsv"),
        _rglob_first(project_root, "dna/gdc-vep.tsv"),
        _rglob_first(project_root, "dna/gdc-dna.tsv"),
        _rglob_first(project_root, "gdc-vep.tsv"),
        _rglob_first(project_root, "gdc-dna.tsv"),
    ],
}

tables = {}
for mod, plist in candidates_map.items():
    pick = _first_existing(plist)
    if pick:
        try:
            tables[mod] = _read_gdc(pick)
            logger.info("Loaded %s table: %s (rows=%d, cols=%d)",
                        mod, pick, *tables[mod].shape)
        except Exception as e:
            logger.warning("Failed reading %s: %s", pick, e)
    else:
        logger.info("No %s GDC table found", mod)

    # Collect union of normalized Case IDs across available tables
    case_ids = set()
    for df in tables.values():
        case_col = "Case ID" if "Case ID" in df.columns else df.columns[0]
        for v in df[case_col].astype(str).tolist():
            v2 = _norm_case_id(v)
            if v2 and v2.lower() != "nan":
                case_ids.add(v2)

    rows = []
    for cid in sorted(case_ids):
        row = {"Case-ID": cid, "DNA": None, "RNA": None, "CH3": None, "CN": None, "Protein": None}

        def pick_path(df, modality):
            if df is None or df.empty:
                return None
            # Normalize Case IDs in a temp column and select
            if "Case ID" in df.columns:
                tmp = df.copy()
                tmp["__CID__"] = tmp["Case ID"].astype(str).map(_norm_case_id)
                df_cid = tmp[tmp["__CID__"] == cid]
            else:
                df_cid = pd.DataFrame()

            if df_cid.empty:
                return None
            pick = df_cid
            if "Sample Type" in df_cid.columns:
                pref = df_cid[df_cid["Sample Type"].astype(str).str.contains(prefer_sample_type, na=False)]
                if not pref.empty:
                    pick = pref
            r = pick.iloc[0]
            fid, fname = r["File ID"], r["File Name"]
            return str(project_root / modality / str(fid) / str(fname))

        row["RNA"] = pick_path(tables.get("rna"), "rna")
        row["CH3"] = pick_path(tables.get("ch3"), "ch3")

        # CN: accept source under cn/ or cnv/; prefer cnv path if it exists, else cn
        cn_df = tables.get("cn")
        if cn_df is not None:
            path_cn = pick_path(cn_df, "cnv")
            if path_cn and not Path(path_cn).exists():
                alt = pick_path(cn_df, "cn")
                path_cn = alt if alt else path_cn
            row["CN"] = path_cn

        # DNA from VEP (preferred) or dna/
        dna_df = tables.get("dna")
        if dna_df is not None:
            path_dna = pick_path(dna_df, "vep")
            if path_dna and not Path(path_dna).exists():
                alt = pick_path(dna_df, "dna")
                path_dna = alt if alt else path_dna
            row["DNA"] = path_dna

        # Protein by filename match
        prot_root = project_root / "protein"
        if prot_root.exists():
            exact = prot_root / f"{cid}.csv"
            if exact.exists():
                row["Protein"] = str(exact)
            else:
                for p in prot_root.glob("*.csv"):
                    if cid in p.stem:
                        row["Protein"] = str(p)
                        break

        rows.append(row)

    out_df = pd.DataFrame(rows, columns=["Case-ID", "DNA", "RNA", "CH3", "CN", "Protein"])
    _write_tsv(out_df, out_path)
    logger.info("Manifest written to %s", out_path)
    return out_path

# ---------- Emit per-modality GDC-like manifests from a unified main manifest ----------

def emit_gdc_like_from_main(main_manifest: Path, project_root: Path, out_dir: Path) -> dict:
    """Create rna_manifest.tsv, ch3_manifest.tsv, cn_manifest.tsv, dna_manifest.tsv under out_dir."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    df = _read_any(Path(main_manifest))
    colmap = {_norm_col(c): c for c in df.columns}

    # Required id column (hyphen-normalized)
    id_col = None
    for key in ("case-id", "caseid"):
        if key in colmap:
            id_col = colmap[key]
            break
    if not id_col:
        raise SystemExit("Main manifest must include a 'Case-ID' (or equivalent) column")

    def as_rel_parts(path_str: str, modality: str):
        pth = Path(path_str)
        if not pth.is_absolute():
            pth = (Path(project_root) / modality / pth).resolve()
        mod_root = (Path(project_root) / modality).resolve()
        try:
            rel = pth.relative_to(mod_root)
            sub = str(rel.parent).replace("\\", "/")
            if sub == ".":
                sub = "."
            fname = rel.name
            return sub, fname
        except Exception:
            return ".", pth.name

    def write_mod(modality: str, col_key: str, outf: str, data_category: str, data_type: str):
        key = _norm_col(col_key)
        if key not in colmap:
            return None
        path_col = colmap[key]
        rows = []
        for _, row in df.iterrows():
            case_id = str(row[id_col])
            path = row[path_col]
            if pd.isna(path) or str(path).strip() == "":
                continue
            sub, fname = as_rel_parts(str(path), modality)
            file_id = sub if sub != "." else ""
            fname_out = fname
            rows.append({
                "File ID": file_id,
                "File Name": fname_out,
                "Data Category": data_category,
                "Data Type": data_type,
                "Project ID": "LOCAL",
                "Case ID": case_id,
                "Sample ID": case_id,
                "Sample Type": "Primary Tumor",
            })
        if not rows:
            return None
        df_out = pd.DataFrame(rows, columns=[
            "File ID","File Name","Data Category","Data Type","Project ID","Case ID","Sample ID","Sample Type"
        ])
        _write_tsv(df_out, out_dir / outf)
        return out_dir / outf

    wrote = {}
    wrote["rna"] = write_mod("rna", "rna", "rna_manifest.tsv",
                             "Transcriptome Profiling", "Gene Expression Quantification")
    wrote["ch3"] = write_mod("ch3", "ch3", "ch3_manifest.tsv",
                             "DNA Methylation", "Methylation Beta Value")
    wrote["cn"]  = write_mod("cnv", "cn",  "cn_manifest.tsv",
                             "Copy Number Variation", "Gene Level Copy Number")
    # DNA: prefer VEP; fallback to dna/
    dna_written = write_mod("vep", "dna", "dna_manifest.tsv",
                            "Simple Nucleotide Variation", "Annotated Somatic Variants")
    if dna_written is None:
        dna_written = write_mod("dna", "dna", "dna_manifest.tsv",
                                "Simple Nucleotide Variation", "Annotated Somatic Variants")
    wrote["dna"] = dna_written
    return wrote

# ---------- Legacy: emit single-type manifest from one GDC TSV ----------

def emit_single_type_manifest(gdc_manifest: Path, out_dir: Path, type_name: str,
                              folder_id_col: int, file_name_col: int, case_id_col: int) -> Path:
    """Given a GDC-like table, produce {type}_manifest.tsv with columns: Case ID | {TYPE}."""
    df = _read_gdc(gdc_manifest)
    for idx in (folder_id_col, file_name_col, case_id_col):
        if idx < 0 or idx >= len(df.columns):
            raise SystemExit("Column index out of range.")
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    folder_series = df.iloc[:, folder_id_col].astype(str)
    fname_series = df.iloc[:, file_name_col].astype(str)
    case_series = df.iloc[:, case_id_col].astype(str).map(_norm_case_id)

    modality_dir = type_name.lower()
    paths = [f"{modality_dir}/{fid}/{fn}" for fid, fn in zip(folder_series, fname_series)]
    out_df = pd.DataFrame({"Case ID": case_series, type_name.upper(): paths})
    out_path = out_dir / f"{type_name.lower()}_manifest.tsv"
    _write_tsv(out_df, out_path)
    return out_path

# ---------- CLI ----------

def parse_args():
    p = argparse.ArgumentParser(description="Build unified and/or per-type manifests for the multiomic pipeline")
    # Mode A: unified main manifest
    p.add_argument("--build-main", action="store_true",
                   help="Build a unified main-manifest.tsv from per-modality GDC TSVs under --project_root")
    p.add_argument("--project_root", help="Project root containing rna/ch3/cn(v)/vep(dna)/protein subfolders")
    p.add_argument("--out_dir", help="Directory to write output manifests")
    p.add_argument("--prefer-sample-type", default="Primary Tumor",
                   help="Preferred Sample Type string when multiple rows per case (default: Primary Tumor)")
    p.add_argument("--log-level", default="INFO",
                   choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                   help="Logging verbosity (default: INFO)")

    # Mode B: emit per-modality manifests from a unified manifest
    p.add_argument("--emit-ref", action="store_true",
                   help="Emit rna/ch3/cn/dna manifests from --main-manifest into --out_dir")
    p.add_argument("--main-manifest", help="Unified TSV/CSV with columns: Case-ID, RNA, CH3, CN, DNA, Protein (paths)")

    # Mode C: legacy single-type manifest emission
    p.add_argument("--gdc-manifest", help="Path to a GDC manifest TSV (for single-type emit)")
    p.add_argument("--type", help="Modality type to emit (rna, ch3, cn, dna)")
    p.add_argument("--folder-id-col", type=int, default=0)
    p.add_argument("--file-name-col", type=int, default=1)
    p.add_argument("--case-id-col", type=int, default=5)
    return p.parse_args()

def main():
    args = parse_args()

    logging.basicConfig(
       level=getattr(logging, args.log_level.upper(), logging.INFO),
       format="%(asctime)s | %(levelname)s | %(name)s | %(message)s"
    )
    logger.info("Starting make_manifest.")

    # Resolve common paths up-front when available
    out_dir = Path(args.out_dir).expanduser().resolve() if args.out_dir else None
    project_root = Path(args.project_root).expanduser().resolve() if args.project_root else None

    if args.build_main:
        if not args.project_root or not args.out_dir:
            raise SystemExit("--build-main requires --project_root and --out_dir")
        out_dir.mkdir(parents=True, exist_ok=True)
        manifest_path = out_dir / "manifest.tsv"
        rna_manifest_path = out_dir / "rna_manifest.tsv"
        ch3_manifest_path = out_dir / "ch3_manifest.tsv"
        cn_manifest_path = out_dir / "cn_manifest.tsv"

        logger.info("  project_root = %s", project_root)
        logger.info("  out_dir      = %s", out_dir)
        logger.info("  manifest.tsv = %s", manifest_path)

        # Build the canonical main manifest
        built = build_main_manifest(project_root, manifest_path, prefer_sample_type=args.prefer_sample_type)
        logger.info("Wrote unified main manifest: %s", built)

        # Derive modality-specific manifests strictly from main manifest
        n_rna = _write_rna_manifest_from_main(manifest_path, rna_manifest_path)
        if n_rna:
            logger.info("Wrote RNA manifest (%d rows): %s", n_rna, rna_manifest_path)
        n_ch3 = _write_ch3_manifest_from_main(manifest_path, ch3_manifest_path)
        if n_ch3:
            logger.info("Wrote CH3 manifest (%d rows): %s", n_ch3, ch3_manifest_path)
        n_cn = _write_cn_manifest_from_main(manifest_path, cn_manifest_path)
        if n_cn:
            logger.info("Wrote CN manifest (%d rows): %s", n_cn, cn_manifest_path)

    if args.emit_ref:
       if not (args.main_manifest and args.out_dir and args.project_root):
           raise SystemExit("--emit-ref requires --main-manifest, --out_dir, and --project_root")
       out_dir = Path(args.out_dir).expanduser().resolve()
       out_dir.mkdir(parents=True, exist_ok=True)
       wrote = emit_gdc_like_from_main(Path(args.main_manifest), project_root, out_dir)

    # Legacy single-type path
    if args.gdc_manifest and not (args.emit_ref or args.build_main):
        if not (args.out_dir and args.type):
            raise SystemExit("--gdc-manifest mode requires --out_dir and --type")
        outp = emit_single_type_manifest(Path(args.gdc_manifest), out_dir,
                                         args.type, args.folder_id_col, args.file_name_col, args.case_id_col)
        logger.info("Wrote %s", outp)

    if not (args.build_main or args.emit_ref or args.gdc_manifest):
        logger.warning(
            "Nothing to do. Try one of:\n"
            "  --build-main --project_root <ROOT> --out_dir <DIR>\n"
            "  --emit-ref --main-manifest manifest.tsv --project_root <ROOT> --out_ref_dir <DIR>\n"
            "  --gdc-manifest <GDC.tsv> --out_dir <DIR> --type <rna|ch3|cn|dna>"
        )

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        logger.error("Interrupted by user")
        sys.exit(130)
    except SystemExit:
        raise
    except Exception as e:
        logger.exception("Unhandled exception: %s", e)
        sys.exit(1)
