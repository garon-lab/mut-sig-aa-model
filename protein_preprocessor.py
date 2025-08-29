#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PROTEIN PREPROCESSOR

This script creates per-case protein CSVs under 'protein/{Case-ID}.csv' that can be used in multiomic integration. 
Note if the flat directory flag was used in downloading CPTAC data, you can leave the folder section blank on your manifest.
Output csvs will have PSM columns listed per case-id: NP (protein accessions), SEQ (peptide sequence), EV (reliability of match), INT (intensity for the selected channel).

This script processes protein PSM files by:
1. Splitting a sample manifest into channels.
2. Creating index folders for output.
3. Filtering and reorganizing raw data by TMT channel.
4. Rejoining processed files into complete datasets named by case-id.

Dependencies: pandas, numpy, argparse, pathlib

Usage:
python protein_preprocessor.py \
  --in_dir <input directory> \
  --manifest <gdc-manifest> \
  --out_dir <output directory> \
  --channel <channel list (e.g. [all|126|127N|127C|128N|128C|129N|129C|130N|130C|131])> \
  --step <step (e.g., [all|split|prep|join])> \
  --jobs 8

Arguments:
(Required)
   --manifest    GDC-like TSV/CSV with at least Case ID, Channel, Directory, File Name
   --folder      Input directory that contains raw data, format <folder>/P{##}/<File>.f{##}.psm>
   --out_dir     Output directory that will contain <protein/<Case-ID>.csv that can be used in multiomic integration
   --channel     Channel list (should mirror what is found in manifest)
   --step        Step(s) to run: split, prep, join, all 
                 *split writes per-channel manifests to <out_dir>/channel_<CH>.txt
                 *prep reads psm parts and filters by --channel, writes to <out_dir>/<case-id>/part-XX.csv
                 *join concatenates part-*.csv for each case-id into <out_dir>/<case-id>.csv
                 *all runs split, then prep, then join

General:
   --jobs        Controls parallel execution, if not provided, script uses min(8, CPU count)
   --cleanup     Deletes intermediate files
   --dry-run     Shows what files would be deleted by cleanup
   --skip-run    Lets cleanup run without re-running command

Troubleshooting
1. If join says "No parts found", ensure prep ran and that FilefXX.psm were found for each ID.
2. IF you see "Unrecognized channel", confirm the Channel matches one of the listed values.
3. For I/O-heavy runs on fast storage, slightly increasing --jobs over core count can help; reduce if you observe disk contention or high memory.
"""

import argparse
import logging
import os
import re
import shutil
import time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

# Channel → column index in PSMs
CHANNEL_MAP = {
    '126': 22, '127N': 23, '127C': 24, '128N': 25, '128C': 26,
    '129N': 27, '129C': 28, '130N': 29, '130C': 30, '131': 31
}

CASE_ID_RE = re.compile(r"^(?:C3[NL]-\d{5}|[A-Z0-9]{3,}[-_]?[A-Z0-9]{2,})$")

ALLOW_FILE = re.compile(r"part-\d{1,3}\.csv$", re.IGNORECASE)


def parse_channels(arg):
    """
    Accepts 'all', a single channel, comma/space separated string,
    or a list/tuple of strings. Returns (selected, bad, raw_parts)
      - selected: None for 'all', or a set of valid channels
      - bad: list of unrecognized tokens
      - raw_parts: the tokens parsed from input
    """
    def tokenize(x):
        if isinstance(x, (list, tuple, set)):
            tokens = []
            for a in x:
                tokens += re.split(r'[,\s]+', str(a))
            return tokens
        return re.split(r'[,\s]+', str(x))

    parts = [t.strip() for t in tokenize(arg) if str(t).strip()]
    # early: if any token == 'all'
    if any(p.lower() == 'all' for p in parts):
        return None, [p for p in parts if p.lower() not in ('all',) and p not in CHANNEL_MAP], parts

    good = [p for p in parts if p in CHANNEL_MAP]
    bad  = [p for p in parts if p and p not in CHANNEL_MAP]
    return (set(good), bad, parts)

def normalize_plex_name(name: str) -> str:
    if not isinstance(name, str):
        return name
    s = name.strip()
    m = re.fullmatch(r'P(\d{1,2})', s, flags=re.IGNORECASE)
    if m:
        return f"P{int(m.group(1)):02d}"
    return s

def _glob_psm_parts(run_dir: Path):
      # Matches ..._fXX.psm
    return sorted(run_dir.glob("*_f[0-9][0-9].psm"))

def _all_case_outputs_exist(case_ids, out_dir: Path):
    """
    Define what 'done' means for a case. Adjust this function if your final output file names differ.
    For example, if you write combined outputs to out_dir / f"{case}.csv" or ".psm":
    """
    ok = True
    missing = []
    for cid in case_ids:
        # EDIT HERE if your expected output differs:
        expected = out_dir / f"{cid}.csv"
        if not expected.exists():
            ok = False
            missing.append((cid, expected))
    return ok, missing

def _is_case_dir(p: Path) -> bool:
    return p.is_dir() and CASE_ID_RE.match(p.name) is not None

def _list_rel_files(root: Path):
    return [Path(dp, fn).relative_to(root) for dp,_,fns in os.walk(root) for fn in fns]

def removable_case_dir(d: Path) -> bool:
    """True if dir is empty or has only part-*.csv (ignoring dotfiles)."""
    for p in d.rglob("*"):
        if not p.is_file():
            continue
        name = p.name
        if name.startswith("."):       
            continue
        if ALLOW_FILE.fullmatch(name):    
            continue
        return False                      
    return True

# ---------- Flexible manifest + fallback path helpers ----------

def read_manifest_flexible(path: str | Path) -> pd.DataFrame:
    """
    Read a protein manifest, accepting multiple header variants and delimiters.
    Canonical columns: 'Folder' (optional), 'File', 'Case-ID', 'Channel'.
    - Accepts 'ID' or 'Case-ID'
    - Accepts 'File' or 'Filename'
    - Accepts tabs or commas
    - Strips whitespace
    """
    path = Path(path)
    df = pd.read_csv(path, sep=None, engine="python")  # auto-detect delimiter

    def norm(s: str) -> str:
        return re.sub(r'[\s_\-]+', '-', str(s).strip().lower())

    normalized = {c: norm(c) for c in df.columns}
    df = df.rename(columns=normalized)

    # Map synonyms -> canonical
    synonym_map = {
        'id': 'Case-ID', 'case-id': 'Case-ID', 'caseid': 'Case-ID',
        'file': 'File', 'filename': 'File',
        'folder': 'Folder',
        'channel': 'Channel'
    }

    # Coalesce duplicates like case-id, case-id.1, ...
    def coalesce_dupes(base: str) -> None:
        cols = [c for c in df.columns if c == base or c.startswith(base + ".")]
        if not cols:
            return
        if len(cols) == 1:
            tgt = synonym_map.get(cols[0], None)
            if tgt:
                df.rename(columns={cols[0]: tgt}, inplace=True)
            return
        combined = df[cols].bfill(axis=1).iloc[:, 0]
        tgt = synonym_map.get(base, base)
        df[tgt] = combined
        for c in cols:
            if c != tgt:
                df.drop(columns=c, inplace=True, errors="ignore")

    for base in ("case-id", "file", "folder", "channel"):
        coalesce_dupes(base)

    # Final rename pass
    for c in list(df.columns):
        if c in synonym_map and synonym_map[c] not in df.columns:
            df.rename(columns={c: synonym_map[c]}, inplace=True)

    # Strip whitespace
    for col in ("Folder", "File", "Case-ID", "Channel"):
        if col in df.columns:
            df[col] = df[col].astype(str).str.strip()

    # Validate required core columns ('Folder' is optional)
    required_core = ["File", "Case-ID", "Channel"]
    missing = [c for c in required_core if c not in df.columns]
    if missing:
        raise ValueError(f"Manifest missing required column(s): {missing}. Got: {list(df.columns)}")

    keep = [c for c in ("Folder", "File", "Case-ID", "Channel") if c in df.columns]
    return df[keep].copy()

# Alias for the rest of the code
read_manifest = read_manifest_flexible

def _read_psm(ps_path: Path) -> pd.DataFrame:
    # PSMs are typically TSV
    try:
        return pd.read_table(ps_path)
    except Exception:
        return pd.read_csv(ps_path, sep=None, engine='python')

def _write_part(df_part: pd.DataFrame, path: Path) -> None:
    df_part.to_csv(path, sep='\t', index=False)

def _read_part(path: Path) -> pd.DataFrame:
    # Robust read for join
    try:
        df = pd.read_csv(path, sep='\t', engine='python', dtype=str)
        if df.shape[1] >= 4:
            return df
    except Exception:
        pass
    try:
        return pd.read_csv(path, sep=None, engine='python', dtype=str, on_bad_lines='skip')
    except TypeError:
        return pd.read_csv(path, sep=None, engine='python', dtype=str)

def index_dirs(out_dir: str | Path, manifest_path: str | Path):
    logger.info("Creating index directories...")
    df = read_manifest(manifest_path)
    for case_id in df['Case-ID']:
        if not case_id or str(case_id).lower() in {'nan', '<na>'}:
            continue
        Path(out_dir, case_id).mkdir(parents=True, exist_ok=True)

def _process_part_task(task):
    """
    Try folder path first; if missing, search anywhere under root by filename.
    """
    folder_root, folder_name, file_prefix, part_num, ch_key, sample_dir, case_id = task
    try:
        folder_root = Path(folder_root)
        folder_n = normalize_plex_name(folder_name) if folder_name else None
        base = f"{file_prefix}f{part_num:02d}.psm"

        ps_file = None
        if folder_n:
            guess = folder_root / folder_n / base
            if guess.exists():
                ps_file = guess

        if ps_file is None:
            cands = list(folder_root.rglob(base))
            if cands:
                ps_file = cands[0]
            else:
                return ("missing", case_id, part_num, str((folder_root / (folder_n or '') / base)))

        # ---- read once ----
        try:
            df1 = _read_psm(ps_file)
        except Exception as e:
            return ("read_error", case_id, part_num, f"{ps_file}: {e}")

        if df1 is None or df1.empty:
            return ("empty", case_id, part_num, str(ps_file))
        df1 = df1.copy()

        # header normalization helper
        def _norm(s: str) -> str:
            return re.sub(r'[\s_]+', '', str(s).lower())

        # Helper: find a column by candidate names (case/space insensitive)
        def find_col(candidates):
            cols_norm = {_norm(c): c for c in df1.columns}
            for cand in candidates:
                key = _norm(cand)
                if key in cols_norm:
                    return cols_norm[key]
            # partial contains match fallback
            for c in df1.columns:
                if any(_norm(k) in _norm(c) for k in candidates):
                    return c
            return None

        # 1) Protein accession column
        acc_col = find_col([
            "Protein Accession","Master Protein Accessions","Accessions","Protein","Accession","ProteinAccession"
        ])
        if acc_col is None:
            # cautious hard-index fallback
            if df1.shape[1] > 13:
                acc_col = df1.columns[13]
            else:
                return ("read_error", case_id, part_num, f"{ps_file}: could not locate accession column")

        # 2) Peptide sequence column
        seq_col = find_col(["Sequence","Peptide Sequence","Peptide"])
        if seq_col is None:
            seq_col = df1.columns[11] if df1.shape[1] > 11 else None

        # 3) Score / EV column (optional)
        ev_col = find_col(["XCorr","Score","Expectation","EValue","EV"])
        if ev_col is None and df1.shape[1] > 16:
            ev_col = df1.columns[16]

        # 4) Gentle NP filter; if empty, keep all
        s_acc_all = df1[acc_col].astype(str)
        df_np = df1[s_acc_all.str.startswith("NP")]
        if df_np.empty:
            df_np = df1
            s_acc = s_acc_all
        else:
            s_acc = df_np[acc_col].astype(str)

        # 5) Parse columns
        NP = s_acc.str.partition('(')[0]
        if seq_col and seq_col in df_np.columns:
            SEQ = (df_np[seq_col].astype(str)
                     .str.split('+', n=1, regex=False).str.get(1)
                     .fillna('')
                     .str.strip()
                     .str[7:])
        else:
            SEQ = ""

        EV = df_np[ev_col] if (ev_col and ev_col in df_np.columns) else ""

        # 6) Intensity column by channel name
        def find_intensity_col(ch_key: str) -> str | None:
            targets = {
                "126":  ["126","tmt126"],
                "127N": ["127n","127 n","tmt127n"],
                "127C": ["127c","127 c","tmt127c"],
                "128N": ["128n","tmt128n"],
                "128C": ["128c","tmt128c"],
                "129N": ["129n","tmt129n"],
                "129C": ["129c","tmt129c"],
                "130N": ["130n","tmt130n"],
                "130C": ["130c","tmt130c"],
                "131":  ["131","131n","131c","tmt131"],
            }[ch_key]
            for c in df_np.columns:
                lc = _norm(c)
                if any(tok in lc for tok in targets):
                    return c
            return None

        int_col = find_intensity_col(ch_key)
        if not int_col:
            # last resort: index-based fallback if shape allows
            try:
                idx = CHANNEL_MAP[ch_key]
                if idx < df_np.shape[1]:
                    int_col = df_np.columns[idx]
            except Exception:
                pass
        if not int_col or int_col not in df_np.columns:
            return ("channel_index_error", case_id, part_num,
                    f"{ps_file}: could not locate intensity for {ch_key}")

        OV = df_np[int_col]

        # 7) Build output
        df_out = pd.DataFrame({"NP": NP.values, "SEQ": SEQ if isinstance(SEQ, pd.Series) else ["" ]*len(NP),
                               "EV": EV if isinstance(EV, pd.Series) else ["" ]*len(NP),
                               "INT": OV.values})
        if df_out.empty:
            return ("empty", case_id, part_num, str(ps_file))

        part_file = Path(sample_dir) / f"part-{part_num}.csv"
        _write_part(df_out, part_file)
        return ("ok", case_id, part_num, str(part_file))

    except Exception as e:
        return ("unexpected_error", case_id, part_num, str(e))



def prep_data(folder, manifest_path, out_dir, channel_flag, jobs=None):
    """
    channel_flag:
      - None  => 'all'
      - set[str] of specific channels (subset of CHANNEL_MAP)
    """
    logger.info("Filtering and organizing PSMs by channel...")
    df = read_manifest(manifest_path)

    # Log selection
    if channel_flag is None:
        logger.info(f"Channel selection: all ({list(CHANNEL_MAP.keys())})")
    else:
        sel_sorted = sorted(channel_flag, key=lambda x: list(CHANNEL_MAP).index(x))
        logger.info(f"Channel selection ({len(sel_sorted)}): {sel_sorted}")

    tasks = []
    kept_rows = skipped_unknown = skipped_not_selected = 0

    for _, row in df.iterrows():
        case_id = row['Case-ID']
        CH = str(row['Channel']).strip()
        FOLDER = (row['Folder'] if 'Folder' in row and pd.notna(row['Folder']) else None)
        FILE = str(row['File']).strip()

        if not case_id or str(case_id).lower() in {'nan', '<na>'}:
            continue
        if CH not in CHANNEL_MAP:
            skipped_unknown += 1
            continue
        if channel_flag is not None and CH not in channel_flag:
            skipped_not_selected += 1
            continue

        kept_rows += 1
        sample_dir = Path(out_dir) / case_id
        sample_dir.mkdir(parents=True, exist_ok=True)
        for j in range(1, 26):
            tasks.append((folder, FOLDER, FILE, j, CH, str(sample_dir), case_id))

    logger.info(
        f"Manifest rows → kept: {kept_rows}, "
        f"skipped_unknown_channel: {skipped_unknown}, "
        f"skipped_not_selected: {skipped_not_selected}"
    )

    if not tasks:
        logger.error("No valid channels found after parsing selection.")
        return

    workers = jobs if (jobs and jobs > 0) else min(8, (os.cpu_count() or 4))
    logger.info(f"Launching parallel processing with {workers} workers for {len(tasks)} parts...")
    status_counts = {}
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futures = [ex.submit(_process_part_task, t) for t in tasks]
        for fut in as_completed(futures):
            status, case_id, part_num, info = fut.result()
            status_counts[status] = status_counts.get(status, 0) + 1
            if status in ("missing", "read_error", "channel_index_error", "unexpected_error", "empty"):
                logger.warning(f"[{status}] {case_id} part {part_num}: {info}")

    logger.info("Parallel prep complete. " + ", ".join(f"{k}={v}" for k, v in status_counts.items()))

def join_parts(parts_root, manifest_path, out_dir):
    logger.info("Joining processed parts...")
    df = read_manifest(manifest_path)
    ids = sorted(set(df['Case-ID']))
    for case_id in ids:
        if not case_id or str(case_id).lower() in {'nan', '<na>'}:
            continue
        part_dir = Path(parts_root) / case_id
        parts = [part_dir / f"part-{j}.csv" for j in range(1, 26)]
        dfs = []
        for p in parts:
            if p.exists():
                try:
                    dfs.append(_read_part(p))
                except Exception as e:
                    logger.warning(f"Failed to read part {p}: {e}")
        if dfs:
            final_path = Path(out_dir) / f"{case_id}.csv"
            pd.concat(dfs, ignore_index=True).to_csv(final_path, index=False)
        else:
            logger.warning(f"No parts found for {case_id}")

def split_channels(manifest_path, out_dir):
    logger.info("Splitting manifest into per-channel manifests...")
    df = read_manifest(manifest_path)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # stable sort: known channels first in CHANNEL_MAP order, unknowns at end
    channels = sorted(
        df['Channel'].dropna().astype(str).unique(),
        key=lambda x: list(CHANNEL_MAP).index(x) if x in CHANNEL_MAP else 999
    )

    for ch in channels:
        sub = df[df['Channel'].astype(str) == ch]
        path = out_dir / f"channel_{ch}.txt"
        sub.to_csv(path, sep='\t', index=False)
        logger.info(f"Wrote {len(sub)} rows -> {path}")

def run_cleanup(out_dir: str, dry_run: bool, require_final_csv: bool = False) -> None:
    root = Path(out_dir).expanduser().resolve()
    if not root.exists():
        print(f"[cleanup] out_dir not found: {root}")
        return
    if str(root) == "/":
        raise SystemExit("Refusing to operate on '/'. Pick a specific out_dir.")

    candidates = [d for d in root.iterdir() if d.is_dir() and CASE_ID_RE.match(d.name)]
    print(f"[cleanup] Scanning {len(candidates)} case folders under {root}")

    removed = 0
    for d in sorted(candidates, key=lambda p: p.name):
        if not removable_case_dir(d):
            continue
        if require_final_csv and not (root / f"{d.name}.csv").exists():
            print(f"  skip: {d} (no final {d.name}.csv)")
            continue
        if dry_run:
            print(f"  DRY-RUN would delete: {d}")
        else:
            shutil.rmtree(d)
            print(f"  deleted: {d}")
            removed += 1

    print(f"[cleanup] {'Would delete' if dry_run else 'Deleted'} {removed} case folders.")

def parse_args():
    p = argparse.ArgumentParser()
    # Backward compatibility: allow both --in_dir and --folder
    p.add_argument('--in_dir', dest="in_dir", help="Input directory containing raw data (e.g., .../P##/FilefXX.psm)")
    p.add_argument('--folder', dest="in_dir", help=argparse.SUPPRESS)  # alias for backward compatibility
    p.add_argument('--manifest',
                   help="TSV/CSV with columns: [Folder?] File, Case-ID, Channel")
    p.add_argument('--out_dir', required=True,
                   help="Output directory")
    p.add_argument(
        '--channel',
        help=("Channel(s): 'all', a single (e.g. 130N), "
              "comma-separated (127N,128N,129N), or space-separated (127N 128N 129N).")
    )
    p.add_argument('--step', default="all", help="all | split | prep | join (default: all)")
    p.add_argument('--jobs', type=int, default=None, help="Number of parallel workers (default: min(8,CPU count))")
    p.add_argument("--cleanup", action="store_true",
               help="Clean case-ID folders under out_dir")
    p.add_argument("--dry-run", action="store_true",
               help="Show what would be removed (no deletions)")
    p.add_argument('--skip-run', action='store_true',
               help='Skip split/prep/join steps (useful with --cleanup).')

    return p.parse_args()


def main():
    args = parse_args()
    
    # --- CLEANUP-ONLY SHORT-CIRCUIT ---
    if args.cleanup and args.skip_run:
        run_cleanup(args.out_dir, dry_run=args.dry_run, require_final_csv=True)
        return
      
    if not args.in_dir or not args.manifest or not args.channel:
        raise SystemExit("--in_dir, --manifest, and --channel are required unless you use --cleanup --skip-run.")

    # Normalize channels
    selected, bad, raw = parse_channels(args.channel)
    logger.info(f"--channel raw: {repr(args.channel)}  parsed: {raw}")
    if bad:
        logger.warning(f"Ignoring unrecognized channels: {bad}")
    channel_flag = selected  # None means 'all'

    # Use in_dir (required via either --in_dir or --folder)
    if not args.in_dir:
        raise SystemExit("Error: you must supply --in_dir (or legacy --folder)")
    in_dir = args.in_dir
    manifest = args.manifest
    out = args.out_dir
    jobs = args.jobs
    step = (args.step or "all").lower()

    Path(out).mkdir(parents=True, exist_ok=True)

    if step in ('all', 'split'):
        split_channels(manifest, out)
    if step in ('all', 'prep'):
        index_dirs(out, manifest)
        prep_data(in_dir, manifest, out, channel_flag, jobs=jobs)
    if step in ('all', 'join'):
        join_parts(out, manifest, out)
    
    if args.cleanup:
        run_cleanup(args.out_dir, dry_run=args.dry_run, require_final_csv=True)


if __name__ == '__main__':
    main()

