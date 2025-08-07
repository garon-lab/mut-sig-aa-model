#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Pipeline for processing VEP-annotated Mutect variant files and generating SNP/SNV summaries, mutational signature, and amino acid substitution matrices.

Dependencies: pandas, gunzip

Usage example: 
    python data_prep_summary.py \
        --folder <folder with mutect files> \
        --manifest <sample_ID manifest> \
        --out <output directory> \
        --mode <1 or 2> \
        [--step all|summary|signatures|matrices]

Arguments:
    --folder   Directory containing the input variant files
    --manifest Manifest file with IDs or gdc sample metadata
    --out      Output directory for results
    --mode     Manifest mode 
                 1 = uses csv/tsv where IDs are in first column (auto-detected)
                 2 = uses gdc manifest-style tsv, uses Case-ID from column 6
    --step     Pipeline step to run 
                 all (default), summary, signatures, or matrices
    
Creates directory:
    manifest.txt           - Simplified manifest file
    summary.csv            - Counts of SNP/SNV per sample
    snp-signature.csv      - SNP signature profile
    snv-signature.csv      - SNV signature profile
    prep/                  - Preprocessed VCF text files (sampleID-mutect.txt)
      snp/
        ├── sampleID-snp.csv
        └── matrices/
              └── sampleID.csv
      snv/
        ├── sampleID-snv.csv
        └── matrices/
              └── sampleID.csv


"""

import pandas as pd
import subprocess
from pathlib import Path
import argparse
import logging
import csv

logging.basicConfig(level=logging.INFO)
AA_LIST = list("ACDEFGHIKLMNPQRSTVWY")

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--folder', required=True, help='Path to folder with input VEP files')
    parser.add_argument('--manifest', required=True, help='Path to manifest file')
    parser.add_argument('--out', required=True, help='Output directory')
    parser.add_argument('--mode', type=int, default=1, help='Mode: 1 = list, 2 = gdc tsv')
    parser.add_argument('--step', default='all', help='Optional step control')
    return parser.parse_args()


def create_manifest(manifest_file: Path, out_dir: Path, mode: int):
    df = pd.read_table(manifest_file)
    out_file = out_dir / "manifest.txt"
    logging.info("Creating manifest file...")

    if mode == 1:
        # Auto-detect delimiter (CSV or TSV)
        with open(manifest_file, 'r') as f:
            sample = f.read(1024)
            try:
                dialect = csv.Sniffer().sniff(sample, delimiters=[',', '\t'])
                delimiter = dialect.delimiter
            except csv.Error:
                logging.error("Unable to determine delimiter. Please check your manifest format.")
                raise

        df = pd.read_csv(manifest_file, delimiter=delimiter, header=None)
        simplified = df.iloc[:, 0].astype(str).str.strip()  # IDs from first column of CSV
    elif mode == 2:
        simplified = df.iloc[:, 5].str.split(',', expand=True).iloc[:, 0].str.strip()
    else:
        raise ValueError("Invalid mode for manifest processing")

    simplified.to_csv(out_file, index=False, header=False)


def make_directories(simplified: Path, out_dir: Path):
    df = pd.read_csv(simplified, header=None)
    ids = df.iloc[:, 0]
    for i in ids:
        (out_dir / "prep").mkdir(parents=True, exist_ok=True)
        (out_dir / "snp").mkdir(parents=True, exist_ok=True)
        (out_dir / "snp/matrices").mkdir(parents=True, exist_ok=True)
        (out_dir / "snv").mkdir(parents=True, exist_ok=True)
        (out_dir / "snv/matrices").mkdir(parents=True, exist_ok=True)


def unzip_files(folder: Path, manifest_file: Path):
    df = pd.read_table(manifest_file, header=None)
    for idx, row in df.iterrows():
        uuid = row[0]
        filename = row[1]
        case_id = row[5].split(',')[0].strip()

        full_path = folder / str(uuid) / str(filename)
        new_path = folder / f"{case_id}.vcf.gz"

        if full_path.exists():
            full_path.rename(new_path)
            subprocess.run(["gunzip", str(new_path)], check=False)


def preprocess_mutect(folder: Path, manifest_file: Path, out_dir: Path):
    df = pd.read_table(manifest_file, header=None)
    for idx, row in df.iterrows():
        case_id = row[5].split(',')[0].strip()
        infile = folder / f"{case_id}.vcf"
        outfile = out_dir / "prep" / f"{case_id}.txt"
        if infile.exists():
            with open(infile, 'r') as fin, open(outfile, 'w') as fout:
                for line in fin:
                    if not line.startswith("##"):
                        fout.write(line)


def summarize_variants(simplified: Path, out_dir: Path):
    df = pd.read_csv(simplified, header=None)
    ids = df.iloc[:, 0]
    summary_file = out_dir / "summary.csv"

    with summary_file.open('w') as f:
        f.write("ID,SNP,SNV\n")

    for sample_id in ids:
        try:
            filepath = out_dir / "prep" / f"{sample_id}.txt"
            vcf_df = pd.read_table(filepath, sep='\t', header=None)
            vcf_df.columns = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'NORMAL', 'TUMOR']

            snp_count = vcf_df[(vcf_df['FILTER'].str.contains("alt", na=False)) & (vcf_df['INFO'].str.contains("missense", na=False))].shape[0]
            snv_count = vcf_df[(vcf_df['FILTER'].str.contains("PASS", na=False)) & (vcf_df['INFO'].str.contains("missense", na=False))].shape[0]

            with summary_file.open('a') as f:
                f.write(f"{sample_id},{snp_count},{snv_count}\n")
        except Exception as e:
            logging.warning(f"Skipping {sample_id}: {e}")


def write_signatures(prep_dir: Path, simplified: Path, out_dir: Path, label: str):
    sample_ids = pd.read_csv(simplified, header=None).iloc[:, 0]
    out_file = out_dir / f"{label}-signature.csv"

    header = ['ID', 'SUM', 'CTGA', 'CAGT', 'GCCG', 'ATTA', 'AGTC', 'ACTG']
    all_rows = []

    for sample_id in sample_ids:
        file_path = prep_dir / f"{sample_id}.txt"

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


def extract_mutations(prep_dir: Path, out_dir: Path, simplified: Path, mutation_type: str):
    assert mutation_type in ["snp", "snv"]
    df = pd.read_csv(simplified, header=None)
    ids = df.iloc[:, 0]

    for sample_id in ids:
        try:
            infile = prep_dir / f"{sample_id}.txt"
            vcf_df = pd.read_table(infile, sep='\t', header=None)
            vcf_df.columns = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'NORMAL', 'TUMOR']

            if mutation_type == "snp":
                filtered = vcf_df[vcf_df['FILTER'].str.contains("alt", na=False)]
            else:
                filtered = vcf_df[vcf_df['FILTER'].str.contains("PASS", na=False)]

            info = filtered['INFO'].str.split('|', expand=True)

            if 15 not in info.columns:
                logging.warning(f"Sample {sample_id}: INFO field has fewer than 16 fields - skipping")
                continue
                
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


def write_matrices(out_dir: Path, simplified: Path):
    df = pd.read_csv(simplified, header=None)
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
    args = parse_args()
    args.out = Path(args.out)
    args.manifest = Path(args.manifest)
    
    # 🔍 Check if manifest file exists and is valid
    if not args.manifest.exists():
        logging.error(f"Manifest file not found: {args.manifest}")
        return

    try:
        df_test = pd.read_table(args.manifest, header=None)
        if args.mode == 1 and df_test.shape[1] < 1:
            logging.error("Manifest file must have at least 1 column for mode 1.")
            return
        elif args.mode == 2 and df_test.shape[1] < 6:
            logging.error("Manifest file must have at least 6 columns for mode 2.")
           
            return
    except Exception as e:
        logging.error(f"Failed to read manifest file: {e}")
        return
    
    # Sets paths
    args.out.mkdir(parents=True, exist_ok=True)
    simplified = args.out / "manifest.txt"
    
    # Generates manifest
    create_manifest(args.manifest, args.out, args.mode)
    
    # Pipeline steps
    if args.step in ("all",):
        make_directories(simplified, args.out)
        unzip_files(args.folder, args.manifest)
        preprocess_mutect(args.folder, args.manifest, args.out)

    if args.step in ("all", "summary"):
        summarize_variants(simplified, args.out)

    if args.step in ("all", "signatures"):
        write_signatures(args.out / "prep", simplified, args.out, "snp")
        write_signatures(args.out / "prep", simplified, args.out, "snv")

    if args.step in ("all", "matrices"):
        extract_mutations(args.out / "prep", args.out / "snp", simplified, "snp")
        extract_mutations(args.out / "prep", args.out / "snv", simplified, "snv")
        write_matrices(args.out, simplified)

    if args.step not in ("all", "summary", "signatures", "matrices"):
        parser.error("Invalid --step value. Choose from: all, summary, signatures, matrices.")



if __name__ == "__main__":
    main()
