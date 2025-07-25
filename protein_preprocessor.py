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
    python protein_preprocessor.py <folder> <manifest> <out> <channel> --step [all|channels|index|prep|join]
"""

import os
import sys
import pandas as pd
import numpy as np
import argparse
from pathlib import Path

def channels(folder, manifest, out, channel_prefix):
    os.chdir(folder)
    df = pd.read_table(manifest)
    print('Separating files into channels...')
    for ch in ['126','127','128','129','130','131']:
        df[ch] = df.iloc[:,3].str.contains(ch)
        out_path = os.path.join(out, f"{channel_prefix}-Ch{ch}.txt")
        df[df[ch]].to_csv(out_path, index=False)

def index(folder, manifest, out):
    os.chdir(folder)
    df = pd.read_csv(manifest)
    F1, F2, IDs = df.iloc[:,0], df.iloc[:,1], df.iloc[:,2]
    for i in range(len(F1)):
        os.makedirs(os.path.join(out, IDs[i]), exist_ok=True)

def prep(folder, manifest, out, channel):
    os.chdir(folder)
    df = pd.read_csv(manifest)
    F1, F2, IDs = df.iloc[:,0], df.iloc[:,1], df.iloc[:,2]
    print('Resorting PSMs by ID...')
    for i in range(len(IDs)):
        base_path = os.path.join(folder, F1[i], F2[i])
        out_path = os.path.join(out, IDs[i])
        for j in range(1, 26):
            part_file = os.path.join(out_path, f"part-{j}.csv")
            ps_file = f"{base_path}f{str(j).zfill(2)}.psm"
            try:
                df1 = pd.read_table(ps_file)
                df1['S1'] = df1.iloc[:,13].str.strip().str[:2]
                df1 = df1[df1['S1'].str.contains('NP')]
                NP = df1.iloc[:,13].str.split('(', expand=True)[0]
                EV = df1.iloc[:,16]
                PS = df1.iloc[:,11].str.split('+', expand=True)[1].str.strip().str[7:]

                channel_map = {'126': 22, '127N': 23, '127C': 24, '128N': 25, '128C': 26,
                               '129N': 27, '129C': 28, '130N': 29, '130C': 30, '131': 31}
                OV = df1.iloc[:, channel_map[channel]]

                df_out = pd.concat([NP, PS, EV, OV], axis=1)
                df_out.columns = ['NP', 'SEQ', 'EV', 'INT']
                df_out.to_csv(part_file, index=False)
            except FileNotFoundError:
                print(f"Missing file for ID: {IDs[i]}, part: {j}")
                continue

def join(folder, manifest, out):
    df = pd.read_csv(manifest)
    IDs = df.iloc[:,0]
    print('Joining files...')
    for ID in IDs:
        part_dir = os.path.join(folder, ID)
        final_path = os.path.join(out, f"{ID}.csv")
        parts = [f"part-{i}.csv" for i in range(1, 26)]
        try:
            combined = pd.concat([pd.read_csv(os.path.join(part_dir, part)) for part in parts])
            combined.to_csv(final_path, index=False)
        except FileNotFoundError:
            print(f"Missing parts for ID: {ID}")
            continue

def main():
    parser = argparse.ArgumentParser(description="Protein processing script for TMT channels.")
    parser.add_argument("folder", help="Directory with raw input files")
    parser.add_argument("manifest", help="Input manifest file")
    parser.add_argument("out", help="Output directory")
    parser.add_argument("channel", help="TMT channel (e.g., 126, 127N, etc.)")
    parser.add_argument("--step", choices=["all", "channels", "index", "prep", "join"], default="all", help="Pipeline step to run")
    args = parser.parse_args()

    if args.step in ["all", "channels"]:
        channels(args.folder, args.manifest, args.out, args.channel)
    if args.step in ["all", "index"]:
        index(args.folder, args.manifest, args.out)
    if args.step in ["all", "prepro"]:
        prepro(args.folder, args.manifest, args.out, args.channel)
    if args.step in ["all", "join"]:
        join(args.out, args.manifest, args.out)

if __name__ == "__main__":
    main()
