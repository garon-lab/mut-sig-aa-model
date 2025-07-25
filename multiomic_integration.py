#!/usr/bin/env python3
"""
Multi-Omic Integration Pipeline

This pipeline integrates multiple omics layers (variant, RNA, methylation, protein, and copy number)
for each sample listed in a manifest. Only the final protein files (with copy number added) are retained
in the output directory. 

Note protein files must be preprocessed and in the format {Sample-ID}.csv (see protein_preprocessor.py)

Use --help to print usage information.
"""

import os
import sys
import csv
import argparse
import shutil
import pandas as pd
import numpy as np
import warnings

warnings.simplefilter("ignore")

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Multi-Omic Integration Pipeline."
    )
    parser.add_argument('--manifest', required=True,
                        help='Tab-delimited manifest file with sample IDs')
    parser.add_argument('--input_var_dir', required=True,
                        help='Directory containing variant CSV files named {sample-ID}.csv')
    parser.add_argument('--input_rna_dir', required=True,
                        help='Directory with RNA expression files')
    parser.add_argument('--rna_manifest', required=True,
                        help='RNA manifest table linking sample IDs to RNA files')
    parser.add_argument('--input_ch3_dir', required=True,
                        help='Directory with methylation (CH3) files')
    parser.add_argument('--ch3_manifest', required=True,
                        help='CH3 manifest table linking sample IDs to CH3 files')
    parser.add_argument('--input_protein_dir', required=True,
                        help='Directory with protein annotation files')
    parser.add_argument('--input_cn_dir', required=True,
                        help='Directory with copy number (CNV) files')
    parser.add_argument('--cn_manifest', required=True,
                        help='CNV manifest table linking sample IDs to CNV files')
    parser.add_argument('--out', default='output',
                        help='Directory for final output (default: output)')
    parser.add_argument('--skip_rna', action='store_true',
                        help='Skip the RNA integration step')
    parser.add_argument('--skip_ch3', action='store_true',
                        help='Skip the CH3 integration step')
    parser.add_argument('--skip_protein', action='store_true',
                        help='Skip the protein integration step')
    parser.add_argument('--skip_cn', action='store_true',
                        help='Skip the copy number integration step')
    parser.add_argument('--help_pipeline', action='store_true',
                        help='Show pipeline help and exit')
    return parser.parse_args()


def print_help():
    help_text = """
Multi-Omic Integration Pipeline Help

Usage:
    python multiomic_pipeline.py \
        --manifest MANIFEST_FILE \
        --input_var_dir VARIANT_DIR \
        --input_rna_dir RNA_DIR \
        --rna_manifest RNA_MANIFEST \
        --input_ch3_dir CH3_DIR \
        --ch3_manifest CH3_MANIFEST \
        --input_protein_dir PROTEIN_DIR \
        --input_cn_dir CNV_DIR \
        --cn_manifest CNV_MANIFEST \
        --out OUTPUT_DIR [--skip_rna] [--skip_ch3] [--skip_protein] [--skip_cn]

Arguments:
    --manifest            Tab-delimited file listing sample IDs (first column)
    --input_var_dir       Directory of variant CSVs named {sample-ID}.csv
    --input_rna_dir       Directory of RNA expression files
    --rna_manifest        Table linking sample IDs to RNA file paths
    --input_ch3_dir       Directory of methylation (CH3) files
    --ch3_manifest        Table linking sample IDs to CH3 file paths
    --input_protein_dir   Directory of protein annotation files
    --input_cn_dir        Directory of CNV files
    --cn_manifest         Table linking sample IDs to CNV file paths
    --out                 Directory to write final integrated files
    --skip_*              Flags to skip specific integration steps

Outputs:
    Only the final protein files (with SNV/SNP, RNA, CH3, protein, and CNV) are written to:
        OUTPUT_DIR/{sample-ID}.csv
"""
    print(help_text)


def read_sample_ids(manifest):
    """Read sample IDs from the first column of the manifest."""
    df = pd.read_table(manifest, header=None)
    return df.iloc[:, 0].astype(str).tolist()


def safe_read_csv(path):
    """Read CSV or return None if missing."""
    try:
        return pd.read_csv(path)
    except Exception:
        warnings.warn(f"Could not read file: {path}")
        return None


def merge_with_reference(df_main, df_ref, on_left, on_right):
    """Left-merge two DataFrames on specified keys."""
    return pd.merge(df_main, df_ref, how='left', left_on=on_left, right_on=on_right)


def add_rna(args):
    if args.skip_rna:
        print("Skipping RNA integration step")
        return
    sample_ids = read_sample_ids(args.manifest)
    df_ref = pd.read_table(args.rna_manifest)
    df_ref['ID'] = df_ref.iloc[:, 5].str.split(',').str[0]
    df_ref['Subfolder'] = df_ref.iloc[:, 0]
    df_ref['Filename'] = df_ref.iloc[:, 1].str.strip().str[:-3]
    df_ref = df_ref.set_index('ID')

    for sample_id in sample_ids:
        file1 = os.path.join(args.input_var_dir, f"{sample_id}.csv")
        df1 = safe_read_csv(file1)
        if df1 is None:
            continue
        try:
            sub = df_ref.loc[sample_id, 'Subfolder']
            fname = df_ref.loc[sample_id, 'Filename']
        except KeyError:
            warnings.warn(f"No RNA manifest entry for {sample_id}")
            continue
        file2 = os.path.join(args.input_rna_dir, sub, fname)
        if not os.path.exists(file2):
            warnings.warn(f"RNA file missing for {sample_id}: {file2}")
            continue
        df_rna = pd.read_table(file2, sep='\t', names=['Ensembl', 'Count'])
        df_rna['Ensembl'] = df_rna['Ensembl'].str.split('.').str[0]
        out_path = os.path.join("tmp/rna", f"{sample_id}.csv")
        df_merged = merge_with_reference(df1, df_rna, 'ENSGene', 'Ensembl')
        df_merged.to_csv(out_path, index=False)


def add_ch3(args):
    if args.skip_ch3:
        print("Skipping CH3 integration step")
        return
    sample_ids = read_sample_ids(args.manifest)
    df_ref = pd.read_table(args.ch3_manifest)
    df_ref['ID'] = df_ref.iloc[:, 5].str.split(',').str[0]
    df_ref['Subfolder'] = df_ref.iloc[:, 0]
    df_ref['Filename'] = df_ref.iloc[:, 1]
    df_ref = df_ref.set_index('ID')

    for sample_id in sample_ids:
        file1 = os.path.join("tmp/rna", f"{sample_id}.csv")
        df1 = safe_read_csv(file1)
        if df1 is None:
            continue
        try:
            sub = df_ref.loc[sample_id, 'Subfolder']
            fname = df_ref.loc[sample_id, 'Filename']
        except KeyError:
            warnings.warn(f"No CH3 manifest entry for {sample_id}")
            continue
        file2 = os.path.join(args.input_ch3_dir, sub, fname)
        if not os.path.exists(file2):
            warnings.warn(f"CH3 file missing for {sample_id}: {file2}")
            continue
        df_ch3 = pd.read_table(file2, sep='\t', names=['cg', 'beta'])
        df_ch3.rename(columns={'cg': 'IlmnID', 'beta': 'beta_val'}, inplace=True)
        out_path = os.path.join("tmp/ch3", f"{sample_id}.csv")
        df_merged = pd.merge(df1, df_ch3, how='left', left_on='IlmnID', right_on='IlmnID')
        df_merged.to_csv(out_path, index=False)


def add_protein(args):
    if args.skip_protein:
        print("Skipping protein integration step")
        return
    sample_ids = read_sample_ids(args.manifest)
    for sample_id in sample_ids:
        file1 = os.path.join("tmp/ch3", f"{sample_id}.csv")
        df1 = safe_read_csv(file1)
        if df1 is None:
            continue
        file2 = os.path.join(args.input_protein_dir, f"{sample_id}.csv")
        df2 = safe_read_csv(file2)
        if df2 is None:
            continue
        out_path = os.path.join("tmp/protein", f"{sample_id}.csv")
        df_merged = merge_with_reference(df1, df2, 'To', 'NP')
        df_merged.to_csv(out_path, index=False)


def add_copy_number(args):
    if args.skip_cn:
        print("Skipping copy number integration step")
        return
    sample_ids = read_sample_ids(args.manifest)
    df_ref = pd.read_table(args.cn_manifest)
    df_ref['ID'] = df_ref.iloc[:, 5].str.split(',').str[0]
    df_ref['Subfolder'] = df_ref.iloc[:, 0]
    df_ref['Filename'] = df_ref.iloc[:, 1]
    df_ref = df_ref.set_index('ID')

    for sample_id in sample_ids:
        file1 = os.path.join("tmp/protein", f"{sample_id}.csv")
        df1 = safe_read_csv(file1)
        if df1 is None:
            continue
        try:
            sub = df_ref.loc[sample_id, 'Subfolder']
            fname = df_ref.loc[sample_id, 'Filename']
        except KeyError:
            warnings.warn(f"No CNV manifest entry for {sample_id}")
            continue
        file2 = os.path.join(args.input_cn_dir, sub, fname)
        if not os.path.exists(file2):
            warnings.warn(f"CNV file missing for {sample_id}: {file2}")
            continue
        df2 = pd.read_table(file2, sep='\t')
        df2['copy_number'].replace('', np.nan, inplace=True)
        df2.dropna(subset=['copy_number'], inplace=True)
        df2['EnsemblID'] = df2.iloc[:, 0].str.split('.').str[0]
        df_cn = df2[['EnsemblID', df2.columns[5]]].rename(columns={df2.columns[5]: 'copy_number'})
        out_path = os.path.join("tmp/protein", f"{sample_id}.csv")
        df_merged = pd.merge(df1, df_cn, how='left', left_on='ENSGene', right_on='EnsemblID')
        df_merged.to_csv(out_path, index=False)


def main():
    args = parse_arguments()
    if args.help_pipeline:
        print_help()
        sys.exit(0)

    # Create tmp directories
    for d in ["tmp/rna", "tmp/ch3", "tmp/protein"]:
        os.makedirs(d, exist_ok=True)

    # Run pipeline steps with skip flags
    add_rna(args)
    add_ch3(args)
    add_protein(args)
    add_copy_number(args)

    # Move final integrated files
    final_dir = args.out
    os.makedirs(final_dir, exist_ok=True)
    for fname in os.listdir("tmp/protein"):
        shutil.move(os.path.join("tmp/protein", fname), os.path.join(final_dir, fname))

    # Clean up temporary files
    shutil.rmtree("tmp")


if __name__ == '__main__':
    main()
