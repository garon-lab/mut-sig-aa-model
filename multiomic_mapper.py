#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Multiomic Mapping Pipeline

Adds gene symbols, names, UniProt IDs, and Illumina methylation probe IDs to VEP-annotated Mutect output.

Dependencies: pandas, numpy, argparse, pathlib, logging

Usage:
    python multiomic_mapper.py \
        --folder <variant dir> \
        --manifest <manifest file> \
        --out <output directory> \
        --ref <reference directory> \
        [--step all|prep|gene|uni|ch3]

Arguments:
    --folder    Directory with variant (SNV/SNP) files
    --manifest  CSV file with sample IDs (first column)
    --out       Output directory
    --ref       Directory with mapping reference files
    --step      Pipeline step to run: all (default), prep, gene, uni, or ch3

## Variant Directory

Note input must be VEP-annotated mutect files with the vcf header removed (see data_prep_summary.py)

### Step Options

- `all`: Full pipeline
- `prep`: Creates directories
- `gene`: Adds gene symbols and names
- 'uni': Includes splitting by chromosome, UniProt and protein ID mapping, and joining
- `ch3`: Add methylation probe matches

## Output Structure

```
out/
├── tmp/	# Deleted after run
│   ├── gene/       # Gene symbol added
│   ├── name/      # Gene name added
│   ├── split/    # Split by chromosome
│   ├── split-uni/  # UniProt added
│   ├── split-uniprot/ # Protein IDs added
│   ├── joined/   # Joined CSVs
├── SAMPLEID.csv  # Final output (example)

"""

import pandas as pd
import numpy as np
import os
import argparse
import shutil
from pathlib import Path
import logging


logging.basicConfig(level=logging.INFO)
CHROM_LIST = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]

def parse_args():
    parser = argparse.ArgumentParser(description="Multiomic Mapping Pipeline")
    parser.add_argument('--folder', required=True,
                        help='Directory with SNV/SNP input files')
    parser.add_argument('--manifest', required=True,
                        help='CSV file with sample IDs (first column)')
    parser.add_argument('--out', required=True,
                        help='Output directory')
    parser.add_argument('--ref', required=True,
                        help='Reference directory with mapping files')
    parser.add_argument('--step', choices=['all','prep','gene','uni','ch3'], default='all',
                        help='Pipeline step to run: all (default), prep, gene, uni, or ch3')
    return parser.parse_args()

def make_dirs(ids, out_dir):
    logging.info("Creating directory structure...")
    for i in ids:
        for sub in ["tmp/gene", "tmp/name", "tmp/split", "tmp/split-uni", "tmp/split-uniprot", "tmp/joined"]:
            Path(out_dir, sub, i).mkdir(parents=True, exist_ok=True)

def add_gene_symbol(ids, folder, out_dir, ref_dir):
    logging.info("Adding gene symbols...")
    gene_map = pd.read_csv(Path(ref_dir, "esng_gene-sym.txt"), sep="\t", names=["Gene", "Ensembl"])
    for i in ids:
        try:
            df = pd.read_csv(Path(folder, f"{i}.csv"))
            merged = df.merge(gene_map, how="left", left_on="ENSGene", right_on="Ensembl").iloc[:, :5]
            merged.to_csv(Path(out_dir, "tmp/gene", i, f"{i}.csv"), index=False)
        except Exception as e:
            logging.warning(f"Skipping {i} in gene symbol step: {e}")

def add_gene_name(ids, out_dir, ref_dir):
    logging.info("Adding gene names...")
    name_map = pd.read_csv(Path(ref_dir, "gene-sym_name.txt"), sep="\t", names=["symbol", "name"])
    for i in ids:
        try:
            df = pd.read_csv(Path(out_dir, "tmp/gene", i, f"{i}.csv"))
            merged = df.merge(name_map, how="left", left_on="Gene", right_on="symbol")
            cols = pd.concat([merged.iloc[:, :5], merged.iloc[:, 6]], axis=1)
            cols.to_csv(Path(out_dir, "tmp/name", i, f"{i}.csv"), index=False)
        except Exception as e:
            logging.warning(f"Skipping {i} in gene name step: {e}")

def split_by_chromosome(ids, out_dir):
    logging.info("Splitting files by chromosome...")
    for i in ids:
        try:
            df = pd.read_csv(Path(out_dir, "tmp/name", i, f"{i}.csv")).iloc[1:]
            for chrom in df["#CHROM"].unique():
                df[df["#CHROM"] == chrom].to_csv(Path(out_dir, "tmp/split", i, f"{chrom}.csv"), index=False)
        except Exception as e:
            logging.warning(f"Skipping {i} in split step: {e}")

def add_uniprot_ids(ids, out_dir, ref_dir):
    logging.info("Adding UniProt IDs...")
    for i in ids:
        for chrom in CHROM_LIST:
            try:
                df = pd.read_csv(Path(out_dir, "tmp/split", i, f"{chrom}.csv"))
                map_file = pd.read_csv(Path(ref_dir, f"uni-ensg{chrom}.txt"), sep="\t", names=["From", "To"])
                merged = df.merge(map_file, how="left", left_on="ENSGene", right_on="To").iloc[:, :7]
                merged.to_csv(Path(out_dir, "tmp/split-uni", i, f"{chrom}.csv"), index=False)
            except Exception:
                continue

def add_protein_ids(ids, out_dir, ref_dir):
    logging.info("Adding protein IDs...")
    for i in ids:
        for chrom in CHROM_LIST:
            try:
                df = pd.read_csv(Path(out_dir, "tmp/split-uni", i, f"{chrom}.csv"))
                map_file = pd.read_csv(Path(ref_dir, "uni-np", f"{chrom}.txt"), sep="\t", names=["From", "To"])
                merged = df.merge(map_file, how="left", on="From")
                merged.to_csv(Path(out_dir, "tmp/split-uniprot", i, f"{chrom}.csv"), index=False)
            except Exception:
                continue

def join_chromosomes(ids, out_dir):
    logging.info("Joining chromosomes...")
    for i in ids:
        joined_path = Path(out_dir, "tmp/joined", f"{i}.csv")
        dfs = []
        for chrom in CHROM_LIST:
            file = Path(out_dir, "tmp/split-uniprot", i, f"{chrom}.csv")
            if file.exists():
                dfs.append(pd.read_csv(file))
        if dfs:
            pd.concat(dfs).to_csv(joined_path, index=False)

def add_methylation(ids, out_dir, ref_dir):
    logging.info("Adding methylation probes...")
    map_df = pd.read_csv(Path(ref_dir, "ch3.csv"))
    map_df = pd.DataFrame({
        "Match": map_df.iloc[:, 0],
        "cg": map_df.iloc[:, 2].str.split(";", expand=True)[0]
    })
    for i in ids:
        try:
            df = pd.read_csv(Path(out_dir, "tmp/joined", f"{i}.csv"))
            merged = df.merge(map_df, how="left", left_on="Gene", right_on="Match").iloc[:, :9]
            merged.to_csv(Path(out_dir, f"{i}.csv"), index=False)
        except Exception as e:
            logging.warning(f"Skipping {i} in methylation step: {e}")


def main():
    args = parse_args()
    folder = Path(args.folder)
    out_dir = Path(args.out)
    ref_dir = Path(args.ref)
    ids = pd.read_csv(args.manifest).iloc[:, 0].astype(str).tolist()

    if args.step in ["all", "prep"]:
        make_dirs(ids, out_dir)

    if args.step in ["all", "gene"]:
        add_gene_symbol(ids, folder, out_dir, ref_dir)
        add_gene_name(ids, out_dir, ref_dir)

    if args.step in ["all", "uni"]:
        split_by_chromosome(ids, out_dir)
        add_uniprot_ids(ids, out_dir, ref_dir)
        add_protein_ids(ids, out_dir, ref_dir)
        join_chromosomes(ids, out_dir)

    if args.step in ["all", "ch3"]:
        add_methylation(ids, out_dir, ref_dir, args.suffix)
        
 # Cleanup
    tmp_path = out_dir / "tmp"
    if tmp_path.exists():
        logging.info("Cleaning up temporary files...")
        shutil.rmtree(tmp_path)


if __name__ == "__main__":
    main()
