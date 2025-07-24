#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Pipeline for processing VEP-annotated Mutect variant files and generating SNP/SNV summaries, mutational signature, and amino acid substitution matrices.

Dependencies: pandas gunzip

Usage example: python aa_1.py <folder> <manifest> <out> <x> <y>

<folder> = path to folder with vep-annotated mutect files
<manifest> = path to GDC manifest or list of IDs
<out> = output directory
<x> = 1 if using gdc manifest ends in extension (e.g., .csv), 2 if multiple entries per field
<y> = 1 if need to create manifest/preprocess, 2 if only creating snv/snp summary files

Creates directory:

output/
├── m.txt
├── summary.csv
├── signature.csv
├── prep/
│   └── sampleID-mutect.txt
├── snp/
│   ├── sampleID-snp.csv
│   └── matrices/
│       └── sampleID.csv
├── snv/
│   ├── sampleID-snv.csv
│   └── matrices/
│       └── sampleID.csv


"""

import pandas as pd
import subprocess
from pathlib import Path
import argparse
import logging

logging.basicConfig(level=logging.INFO)
AA_LIST = list("ACDEFGHIKLMNPQRSTVWY")

def print_help():
    help_text = """
    Amino Acid Pipeline Help

    Usage:
        python aa_1.py <folder> <manifest> <out> <x> <y>

    Arguments:
        <folder>   - Directory containing the input variant files
        <manifest> - Tab-delimited manifest file with sample metadata
        <out>      - Output directory to store results
        <x>        - Manifest mode (1 = use column 2 trimmed; 2 = use column 6 split)
        <y>        - Step flag (1 = do full preprocessing; 0 = skip to summary/matrix steps)

    Outputs:
        - m.txt: Simplified manifest file
        - summary.csv: Counts of SNP/SNV per sample
        - Preprocessed VCF text files in prep/
        - Extracted amino acid changes in snp/ and snv/
        - Amino acid substitution matrices in snp/matrices/ and snv/matrices/

    Requirements:
        - Python 3, pandas
        - gunzip must be available for decompressing files
    """


def create_manifest(folder: Path, manifest_file: Path, out_dir: Path, mode: int):
    df = pd.read_table(manifest_file)
    out_file = out_dir / "m.txt"
    logging.info("Creating manifest file...")

    if mode == 1:
        df = df.iloc[:, 1].str.strip().str[:-3]
    elif mode == 2:
        df = df.iloc[:, 5].str.split(',', expand=True).iloc[:, 0]
    else:
        raise ValueError("Invalid mode for manifest processing")

    df.to_csv(out_file, index=False, header=False)


def make_directories(manifest_file: Path, out_dir: Path):
    df = pd.read_csv(manifest_file, header=None)
    ids = df.iloc[:, 0]
    for i in ids:
        (out_dir / i).mkdir(parents=True, exist_ok=True)
        for subfolder in ["prep", "snp/matrices", "snv/matrices"]:
            (out_dir / i / subfolder).mkdir(parents=True, exist_ok=True)


def unzip_files(folder: Path, manifest_file: Path):
    df = pd.read_table(manifest_file, header=None)
    ids = df.iloc[:, 0]
    names = df.iloc[:, 1]
    for i, name in zip(ids, names):
        filepath = folder / i / name
        if filepath.exists():
            subprocess.run(["gunzip", str(filepath)], check=False)


def preprocess_mutect(folder: Path, manifest_file: Path, out_dir: Path):
    df = pd.read_table(manifest_file, header=None)
    filepaths = df.iloc[:, 0]
    names = df.iloc[:, 1].str.strip().str[:-3]
    sample_ids = df.iloc[:, 5].str.split(',', expand=True).iloc[:, 0]

    for path, name, sid in zip(filepaths, names, sample_ids):
        infile = folder / path / name
        outfile = out_dir / "prep" / f"{sid}-mutect.txt"
        if infile.exists():
            with open(infile, 'r') as fin, open(outfile, 'w') as fout:
                for line in fin:
                    if not line.startswith("##"):
                        fout.write(line)


def summarize_variants(out_dir: Path, manifest_file: Path):
    df = pd.read_table(manifest_file, header=None)
    ids = df.iloc[:, 0]
    summary_file = out_dir / "summary.csv"

    with summary_file.open('w') as f:
        f.write("ID,SNP,SNV\n")

    for sample_id in ids:
        try:
            filepath = out_dir / "prep" / f"{sample_id}-mutect.txt"
            vcf_df = pd.read_table(filepath, sep='\t', header=None)
            vcf_df.columns = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'NORMAL', 'TUMOR']

            snp_count = vcf_df[(vcf_df['FILTER'].str.contains("alt", na=False)) & (vcf_df['INFO'].str.contains("missense", na=False))].shape[0]
            snv_count = vcf_df[(vcf_df['FILTER'].str.contains("PASS", na=False)) & (vcf_df['INFO'].str.contains("missense", na=False))].shape[0]

            with summary_file.open('a') as f:
                f.write(f"{sample_id},{snp_count},{snv_count}\n")
        except Exception as e:
            logging.warning(f"Skipping {sample_id}: {e}")
from pathlib import Path
import pandas as pd
import logging


def write_signatures(mutect_dir: Path, manifest_file: Path, out_dir: Path, label: str):
    sample_ids = pd.read_table(manifest_file, header=None).iloc[:, 0]
    out_file = out_dir / f"{label}-signature.csv"
    out_dir.mkdir(parents=True, exist_ok=True)

    header = ['ID', 'SUM', 'CTGA', 'CAGT', 'GCCG', 'ATTA', 'AGTC', 'ACTG']
    all_rows = []

    for sample_id in sample_ids:
        file_path = mutect_dir / f"{sample_id}-mutect.txt"

        if not file_path.exists():
            logging.warning(f"{sample_id}: file not found at {file_path}")
            continue

        try:
            df = pd.read_table(file_path)

            df['sA'] = df.iloc[:, 3].str.contains('A', na=False)
            df['sC'] = df.iloc[:, 3].str.contains('C', na=False)
            df['sG'] = df.iloc[:, 3].str.contains('G', na=False)
            df['sT'] = df.iloc[:, 3].str.contains('T', na=False)
            df['eA'] = df.iloc[:, 4].str.contains('A', na=False)
            df['eC'] = df.iloc[:, 4].str.contains('C', na=False)
            df['eG'] = df.iloc[:, 4].str.contains('G', na=False)
            df['eT'] = df.iloc[:, 4].str.contains('T', na=False)

            df['NM'] = df.iloc[:, 6].str.strip().str[:3]
            df1 = df[df['NM'].str.contains('alt', na=False)]

            # Define mutation combinations
            mutations = {
                'AC': df1['sA'] & df1['eC'],
                'AG': df1['sA'] & df1['eG'],
                'AT': df1['sA'] & df1['eT'],
                'CA': df1['sC'] & df1['eA'],
                'CG': df1['sC'] & df1['eG'],
                'CT': df1['sC'] & df1['eT'],
                'GA': df1['sG'] & df1['eA'],
                'GC': df1['sG'] & df1['eC'],
                'GT': df1['sG'] & df1['eT'],
                'TA': df1['sT'] & df1['eA'],
                'TC': df1['sT'] & df1['eC'],
                'TG': df1['sT'] & df1['eG'],
            }

            counts = {
                'CTGA': (mutations['CT'] | mutations['GA']).sum(),
                'CAGT': (mutations['CA'] | mutations['GT']).sum(),
                'GCCG': (mutations['CG'] | mutations['GC']).sum(),
                'ATTA': (mutations['AT'] | mutations['TA']).sum(),
                'AGTC': (mutations['AG'] | mutations['TC']).sum(),
                'ACTG': (mutations['AC'] | mutations['TG']).sum(),
            }

            total = sum(counts.values())
            row = [sample_id, total] + [counts[key] for key in ['CTGA', 'CAGT', 'GCCG', 'ATTA', 'AGTC', 'ACTG']]
            all_rows.append(row)

        except Exception as e:
            logging.warning(f"Failed to process {sample_id}: {e}")

    # Save all results to CSV
    df_out = pd.DataFrame(all_rows, columns=header)
    df_out.to_csv(out_file, index=False)


def extract_mutations(prep_dir: Path, out_dir: Path, manifest_file: Path, mutation_type: str):
    assert mutation_type in ["snp", "snv"]
    df = pd.read_table(manifest_file, header=None)
    ids = df.iloc[:, 0]

    for sample_id in ids:
        try:
            infile = prep_dir / f"{sample_id}-mutect.txt"
            vcf_df = pd.read_table(infile, sep='\t', header=None)
            vcf_df.columns = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'NORMAL', 'TUMOR']

            if mutation_type == "snp":
                filtered = vcf_df[vcf_df['FILTER'].str.contains("alt", na=False)]
            else:
                filtered = vcf_df[vcf_df['FILTER'].str.contains("PASS", na=False)]

            info = filtered['INFO'].str.split('|', expand=True)
            aa_df = pd.DataFrame({
                'ST': info[15].str[0],
                'END': info[15].str[-1],
                '#CHROM': filtered['#CHROM'],
                'TUMOR': filtered['TUMOR']
            })
            outfile = out_dir / f"{sample_id}-{mutation_type}.csv"
            aa_df.to_csv(outfile, index=False)
        except Exception as e:
            logging.warning(f"Failed to extract {mutation_type} for {sample_id}: {e}")


def generate_aa_matrix(df):
    matrix = df.groupby(['ST', 'END']).size().unstack(fill_value=0)
    matrix = matrix.reindex(index=AA_LIST + ['*'], columns=AA_LIST + ['*'], fill_value=0)
    return matrix


def write_matrices(out_dir: Path, manifest_file: Path):
    df = pd.read_table(manifest_file, header=None)
    ids = df.iloc[:, 0]
    for sample_id in ids:
        for mtype in ['snp', 'snv']:
            try:
                file = out_dir / mtype / f"{sample_id}-{mtype}.csv"
                if not file.exists():
                    continue
                df_aa = pd.read_csv(file)
                matrix = generate_aa_matrix(df_aa)
                out_path = out_dir / mtype / "matrices" / f"{sample_id}.csv"
                matrix.to_csv(out_path)
            except Exception as e:
                logging.warning(f"Matrix generation failed for {sample_id}-{mtype}: {e}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("folder", type=Path)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("out", type=Path)
    parser.add_argument("x", type=int)
    parser.add_argument("y", type=int)
    parser.add_argument("--help-pipeline", action="store_true", help="Display usage information for the pipeline")
    args = parser.parse_args()

    args.out.mkdir(parents=True, exist_ok=True)

    if args.y == 1:
        create_manifest(args.folder, args.manifest, args.out, args.x)
        make_directories(args.out / "m.txt", args.out)
        unzip_files(args.folder, args.manifest)
        preprocess_mutect(args.folder, args.manifest, args.out)

    summarize_variants(args.out, args.manifest)
    extract_mutations(args.out / "prep", args.out / "snp", args.manifest, "snp")
    extract_mutations(args.out / "prep", args.out / "snv", args.manifest, "snv")
    write_matrices(args.out, args.manifest)


if __name__ == "__main__":
    main()
