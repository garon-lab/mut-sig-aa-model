#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Protein Processing Pipeline

This script processes protein PSM files by:
1. Splitting a sample manifest into channels.
2. Creating index folders for output.
3. Filtering and reorganizing raw data by TMT channel.
4. Rejoining processed files into complete datasets.

Usage:
    python protein_preprocessor.py \
        --folder <input directory> \
        --manifest <manifest file> \
        --out <output directory> \
        --channel <TMT channel (e.g., 126,127N)> \
        [--step all|channels|index|prep|join]

Arguments:
    --folder    Directory with raw input files
    --manifest  Tab-delimited manifest file listing sample metadata
    --out       Output directory to store results
    --channel   TMT channel identifier (e.g., 126, 127N, etc.)
    --step      Pipeline step to run: all (default), channels, index, prep, or join
"""

import os
import sys
import pandas as pd
import numpy as np
import argparse
from pathlib import Path
import logging

logging.basicConfig(level=logging.INFO)


def parse_args():
    parser = argparse.ArgumentParser(description="Protein processing pipeline for TMT channels.")
    parser.add_argument('--folder', required=True,
                        help='Directory containing raw input PSM files')
    parser.add_argument('--manifest', required=True,
                        help='Tab-delimited manifest file with sample metadata')
    parser.add_argument('--out', required=True,
                        help='Output directory to store processed results')
    parser.add_argument('--channel', required=True,
                        help='TMT channel to process (e.g., 126,127N,128C)')
    parser.add_argument('--step', choices=['all', 'channels', 'index', 'prep', 'join'],
                        default='all',
                        help='Pipeline step to run (default: all)')
    return parser.parse_args()


def split_channels(folder, manifest, out, channel_prefix):
    os.chdir(folder)
    df = pd.read_table(manifest)
    logging.info('Splitting manifest into channels...')
    for ch in ['126','127','128','129','130','131']:
        df[ch] = df.iloc[:,3].str.contains(ch)
        out_path = Path(out) / f"{channel_prefix}-Ch{ch}.txt"
        df[df[ch]].to_csv(out_path, index=False)


def index_dirs(folder, manifest, out):
    os.chdir(folder)
    df = pd.read_csv(manifest)
    IDs = df.iloc[:,2]
    logging.info('Creating index directories...')
    for i in IDs:
        Path(out, i).mkdir(parents=True, exist_ok=True)


def prep_data(folder, manifest, out, channel):
    os.chdir(folder)
    df = pd.read_csv(manifest)
    F1, F2, IDs = df.iloc[:,0], df.iloc[:,1], df.iloc[:,2]
    logging.info('Filtering and organizing PSMs by channel...')
    # map channel codes to intensity column
    channel_map = {'126':22,'127N':23,'127C':24,'128N':25,'128C':26,
                   '129N':27,'129C':28,'130N':29,'130C':30,'131':31}
    for i, f1, f2 in zip(IDs, F1, F2):
        sample_dir = Path(out) / i
        for j in range(1, 26):
            ps_file = Path(folder) / f1 / f2 / f"f{str(j).zfill(2)}.psm"
            part_file = sample_dir / f"part-{j}.csv"
            if not ps_file.exists():
                logging.warning(f"Missing PSM: {i} part {j}")
                continue
            df1 = pd.read_table(ps_file)
            df1['S1'] = df1.iloc[:,13].str[:2]
            df1 = df1[df1['S1'] == 'NP']
            NP = df1.iloc[:,13].str.split('(',1).str[0]
            EV = df1.iloc[:,16]
            PS = df1.iloc[:,11].str.split('+',1).str[1].str.strip().str[7:]
            OV = df1.iloc[:, channel_map[channel]]
            df_out = pd.DataFrame({'NP':NP, 'SEQ':PS, 'EV':EV, 'INT':OV})
            df_out.to_csv(part_file, index=False)


def join_parts(folder, manifest, out):
    df = pd.read_csv(manifest)
    IDs = df.iloc[:,0]
    logging.info('Joining processed parts...')
    for ID in IDs:
        part_dir = Path(folder) / ID
        final_path = Path(out) / f"{ID}.csv"
        parts = [part_dir / f"part-{j}.csv" for j in range(1,26)]
        dfs = []
        for p in parts:
            if p.exists():
                dfs.append(pd.read_csv(p))
        if dfs:
            pd.concat(dfs).to_csv(final_path, index=False)
        else:
            logging.warning(f"No parts found for {ID}")


def main():
    args = parse_args()
    folder = args.folder
    manifest = args.manifest
    out = args.out
    channel = args.channel

    if args.step in ['all', 'channels']:
        split_channels(folder, manifest, out, channel)
    if args.step in ['all', 'index']:
        index_dirs(folder, manifest, out)
    if args.step in ['all', 'prep']:
        prep_data(folder, manifest, out, channel)
    if args.step in ['all', 'join']:
        join_parts(out, manifest, out)

if __name__ == "__main__":
    main()
