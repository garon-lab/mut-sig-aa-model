# MUT-SIG-AA-MODEL
Uses publicly available genomic data (e.g., CPTAC, TCGA) to create models of amino acid substitutions from non-synonymous mutations with ability to compare amino acid substitution matrices. 

Order of scripts:
1. Data prep summary
2. Protein (PSM) preprocesser
3. Multiomic mapping
4. Multiomic data integration
5. Analysis

Additional information regarding each script is below and can be accessed with the --help flag.


# DATA PREP SUMMARY

Pipeline for processing VEP-annotated Mutect variant files and generating SNP/SNV summaries, mutational signature, and amino acid substitution matrices.

Dependencies: pandas, gunzip

Usage example: 
python data_prep_summary.py \
        --folder <folder with mutect files> \
        --manifest <sample_ID manifest> \
        --out <output directory> \
        --mode <mode> \
        [--step all|summary|signatures|matrices]

Arguments:
    --folder   Directory containing the input variant files
    --manifest Tab-delimited manifest file with sample metadata
    --out      Output directory for results
    --mode     Manifest mode (1 = use column 2 trimmed (removes .gz from naming); 2 = use column 6 split for multiple samples)
    --step     Pipeline step to run: all (default), summary, signatures, or matrices
    
Creates directory:
    m.txt                  - Simplified manifest file
    summary.csv            - Counts of SNP/SNV per sample
    signature.csv          - Mutational signature files for SNP/SNV
    prep/                  - Preprocessed VCF text files (sampleID-mutect.txt)
    snp/                   - Extracted SNP amino acid changes
        matrices/          - SNP substitution matrices
    snv/                   - Extracted SNV amino acid changes
        matrices/          - SNV substitution matrices


# PROTEIN PREPROCESSOR

This script processes protein PSM files by:
1. Splitting a sample manifest into channels.
2. Creating index folders for output.
3. Filtering and reorganizing raw data by TMT channel.
4. Rejoining processed files into complete datasets named by sample-id.

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


# MULTIOMIC MAPPER

This pipeline augments VEP-annotated Mutect variant files with gene symbols, gene names, UniProt IDs, and Illumina methylation probe IDs ('cg'), splitting by chromosome and recombining to enable the addition of multiomic data to to VEP-annotated Mutect files.

Dependencies: pandas, numpy, argparse, pathlib, logging

Usage example: 
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

Step options:
- all: Full pipeline
- prep: Creates directories
- gene: Adds gene symbols and names
- uni: Includes splitting by chromosome, UniProt and protein ID mapping, and joining
- ch3: Add methylation probe matches


# MULTIOMIC INTEGRATION

This pipeline integrates multiple omics layers (variant, RNA, methylation, protein, and copy number)
for each sample listed in a manifest. Only the final protein files (with copy number added) are retained
in the output directory. Note protein files must be preprocessed and in the format {Sample-ID}.csv (see protein_preprocessor.py)

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
