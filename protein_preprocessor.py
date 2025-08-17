#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PROTEIN PREPROCESSOR

This script creates per-case protein CSVs under 'protein/{Case-ID}.csv' that can be used in multiomic integration. Output csvs will have PSM columns listed per case-id: NP (protein accessions), SEQ (peptide sequence), EV (reliability of match), INT (intensity for the selected channel).

This script processes protein PSM files by:
1. Splitting a sample manifest into channels.
2. Creating index folders for output.
3. Filtering and reorganizing raw data by TMT channel.
4. Rejoining processed files into complete datasets named by case-id.

Dependencies: pandas, numpy, argparse, pathlib

Usage:
python protein_preprocessor.py \
  --folder <input directory> \
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
                 *join concatenantes part-*.csv for each case-id into <out_dir>/<case-id>.csv
                 *all runs split, then prep, then join

General:
   --jobs        Controls parallel execution, if not provided, script uses min(8, CPU count)

Troubleshooting
1. If join says "No parts found", ensure prep ran and that FilefXX.psm were found for each ID.
2. IF you see "Unrecognizaed channel", confirm the Channel matches one of the listed values.
3. For I/O-heavy runs on fast storage, slightly increasing --jobs over core count can help; reduce if you observe disk contention or high memory.
"""
import argparse
import logging
import os
import re
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

import pandas as pd

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Channel → column index in PSMs
CHANNEL_MAP = {'126':22,'127N':23,'127C':24,'128N':25,'128C':26,
               '129N':27,'129C':28,'130N':29,'130C':30,'131':31}

def normalize_plex_name(name: str) -> str:
    if not isinstance(name, str):
        return name
    s = name.strip()
    m = re.fullmatch(r'P(\d{1,2})', s, flags=re.IGNORECASE)
    if m:
        return f"P{int(m.group(1)):02d}"
    return s

def read_manifest(manifest):
    df = pd.read_csv(manifest, sep=None, engine='python', dtype=str, keep_default_na=False)
    cols = [c.strip() for c in df.columns.tolist()]
    df.columns = cols
    need = {'Case-ID','Channel','Folder','File'}
    # Back-compat: if legacy 'ID' present, rename to 'Case-ID'
    if 'Case-ID' not in df.columns and 'ID' in df.columns:
        df = df.rename(columns={'ID':'Case-ID'})
    if not need.issubset(set(cols)):
        if len(cols) >= 4:
            rename_map = {cols[0]:'Case-ID', cols[1]:'Channel', cols[2]:'Folder', cols[3]:'File'}
            df = df.rename(columns=rename_map)
    for c in ['Case-ID','Channel','Folder','File']:
        if c not in df.columns:
            raise ValueError(f"Manifest missing required column '{c}'. Got: {df.columns.tolist()}")
    for c in ['Case-ID','Channel','Folder','File','Type']:
        if c in df.columns:
            df[c] = df[c].astype(str).str.strip()
    return df

def index_dirs(folder, manifest, out_dir):
    logger.info("Creating index directories...")
    df = read_manifest(manifest)
    for case_id in df['Case-ID']:
        if not case_id or str(case_id).lower() in {'nan','<na>'}:
            continue
        Path(out_dir, case_id).mkdir(parents=True, exist_ok=True)

def _read_psm(ps_path: Path) -> pd.DataFrame:
    # PSMs are TSV; be forgiving if not
    try:
        return pd.read_table(ps_path)
    except Exception:
        return pd.read_csv(ps_path, sep=None, engine='python')

def _write_part(df_part: pd.DataFrame, path: Path) -> None:
    # Write TSV to avoid comma/quote issues
    df_part.to_csv(path, sep='\t', index=False)

def _read_part(path: Path) -> pd.DataFrame:
    # Robust read for join: try TSV then sniff; tolerate bad lines
    try:
        df = pd.read_csv(path, sep='\t', engine='python', dtype=str)
        if df.shape[1] >= 4:
            return df
    except Exception:
        pass
    try:
        df = pd.read_csv(path, sep=None, engine='python', dtype=str, on_bad_lines='skip')
        return df
    except TypeError:
        # older pandas: on_bad_lines might not exist
        df = pd.read_csv(path, sep=None, engine='python', dtype=str)
        return df

def _process_part_task(task):
    folder_root, folder_name, file_prefix, part_num, ch_key, sample_dir, case_id = task
    try:
        folder_n = normalize_plex_name(folder_name)
        base = f"{file_prefix}f{part_num:02d}.psm"
        ps_file = Path(folder_root) / folder_n / base
        if not ps_file.exists():
            # Glob fallback
            folder_dir = Path(folder_root) / folder_n
            cands = sorted(folder_dir.glob(f"{file_prefix}*f{part_num:02d}.psm"))
            if cands:
                ps_file = cands[0]
            else:
                return ("missing", case_id, part_num, str(ps_file))

        try:
            df1 = _read_psm(ps_file)
        except Exception as e:
            return ("read_error", case_id, part_num, f"{ps_file}: {e}")

        # Filter and extract
        df1 = df1.copy()
        df1['S1'] = df1.iloc[:,13].astype(str).str[:2]
        df1 = df1[df1['S1'] == 'NP']
        if df1.empty:
            return ("empty", case_id, part_num, str(ps_file))

        NP = df1.iloc[:,13].astype(str).str.partition('(')[0]
        EV = df1.iloc[:,16]
        PS = df1.iloc[:,11].astype(str).str.split('+', n=1, regex=False).str.get(1).fillna('').str.strip().str[7:]
        try:
            OV = df1.iloc[:, CHANNEL_MAP[ch_key]]
        except Exception as e:
            return ("channel_index_error", case_id, part_num, f"{ps_file}: {e}")

        df_out = pd.DataFrame({'NP':NP, 'SEQ':PS, 'EV':EV, 'INT':OV})
        part_file = Path(sample_dir) / f"part-{part_num}.csv"
        _write_part(df_out, part_file)
        return ("ok", case_id, part_num, str(part_file))
    except Exception as e:
        return ("unexpected_error", case_id, part_num, str(e))

def prep_data(folder, manifest, out_dir, channel_flag, jobs=None):
    logger.info("Filtering and organizing PSMs by channel...")
    df = read_manifest(manifest)

    tasks = []
    for _, row in df.iterrows():
        case_id = row['Case-ID']; CH = row['Channel']; FOLDER = row['Folder']; FILE = row['File']
        if not case_id or str(case_id).lower() in {'nan','<na>'}:
            continue
        # honor --channel all vs explicit channel
        if isinstance(channel_flag, str) and channel_flag.lower() == 'all':
            ch_key = CH
        else:
            ch_key = channel_flag
        if CH and ch_key and (isinstance(channel_flag, str) and channel_flag.lower() != 'all') and (CH != ch_key):
            continue
        if ch_key not in CHANNEL_MAP:
            logger.warning(f"[prep] Unrecognized channel '{ch_key}' for Case-ID {case_id}; skipping sample.")
            continue
        sample_dir = Path(out_dir) / case_id
        sample_dir.mkdir(parents=True, exist_ok=True)
        for j in range(1, 26):
            tasks.append((folder, FOLDER, FILE, j, ch_key, str(sample_dir), case_id))

    if not tasks:
        logger.warning("No tasks to run (check manifest or channel selection).")
        return

    workers = jobs if (jobs and jobs > 0) else (os.cpu_count() or 4)
    logger.info(f"Launching parallel processing with {workers} workers for {len(tasks)} parts...")
    status_counts = {}
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futures = [ex.submit(_process_part_task, t) for t in tasks]
        for fut in as_completed(futures):
            status, case_id, part_num, info = fut.result()
            status_counts[status] = status_counts.get(status, 0) + 1
            if status in ("missing","read_error","channel_index_error","unexpected_error","empty"):
                logger.warning(f"[{status}] {case_id} part {part_num}: {info}")

    logger.info("Parallel prep complete. Status summary: " + ", ".join(f"{k}={v}" for k,v in status_counts.items()))

def join_parts(parts_root, manifest, out_dir):
    logger.info("Joining processed parts...")
    df = read_manifest(manifest)
    ids = sorted(set(df['Case-ID']))
    for case_id in ids:
        if not case_id or str(case_id).lower() in {'nan','<na>'}:
            continue
        part_dir = Path(parts_root) / case_id
        parts = [part_dir / f"part-{j}.csv" for j in range(1,26)]
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

def split_channels(folder, manifest, out_dir, channel_flag):
    logger.info("Splitting manifest into channels (optional utility)...")
    df = read_manifest(manifest)
    out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)
    for ch in CHANNEL_MAP.keys():
        path = out_dir / f"channel_{ch}.txt"
        df[df['Channel'] == ch].to_csv(path, index=False)

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('--folder', required=True)
    p.add_argument('--manifest', required=True)
    p.add_argument('--out_dir', required=True)
    p.add_argument('--channel', required=True, help="e.g., 126,127N,... or 'all'")
    p.add_argument('--step', required=True, help="all | split | prep | join")
    p.add_argument('--jobs', type=int, default=None)
    return p.parse_args()

def main():
    args = parse_args()
    folder = args.folder
    manifest = args.manifest
    out = args.out_dir
    step = args.step.lower()
    channel_flag = args.channel
    jobs = args.jobs

    if step in ('all','split'):
        split_channels(folder, manifest, out, channel_flag)
    if step in ('all','prep'):
        index_dirs(folder, manifest, out)
        prep_data(folder, manifest, out, channel_flag, jobs=jobs)
    if step in ('all','join'):
        # IMPORTANT: join_parts expects its first argument to be the folder where parts live (== out)
        join_parts(out, manifest, out)

if __name__ == '__main__':
    main()
