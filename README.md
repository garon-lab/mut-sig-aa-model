# MUT-SIG-AA-MODEL
Uses publicly available genomic data (e.g., CPTAC, TCGA) to create models of amino acid substitutions from non-synonymous mutations with ability to compare amino acid substitution matrices

# DATA PREP SUMMARY

Pipeline for processing VEP-annotated Mutect variant files and generating SNP/SNV summaries, mutational signature, and amino acid substitution matrices.

Dependencies: pandas gunzip

Usage example: python data_prep_summary.py <folder> <manifest> <out> <x> <y>

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


# MULTIOMIC MAPPER

This pipeline augments VEP-annotated Mutect variant files with gene symbols, gene names, UniProt IDs, and Illumina methylation probe IDs ('cg'), splitting by chromosome and recombining for analysis.

## Requirements

- Python 3
- pandas
- numpy

## Usage

bash
python multiomic_mapper.py <folder> <manifest> <out> <ref> <suffix> --step [all|prep|gene|uni|methyl]


### Arguments

- <folder>: Input folder with VEP-annotated SNV/SNP files
- <manifest>: CSV with first column = sample IDs
- <out>: Output directory
- <ref>: Folder with mapping references (e.g., Ensembl-symbol, UniProt-ENSG, cg IDs)
- <suffix>: File suffix (e.g., '-snv.csv')
- --step: Which step to run (default 'all')

### Step Options

- all: Full pipeline
- prep: Creates directories
- gene: Adds gene symbols and names
- uni: Includes splitting by chromosome, UniProt and protein ID mapping, and joining
- ch3: Add methylation probe matches

## Output Structure

out/
├── tmp/	# Deleted after run
│   ├── gene/       # Gene symbol added
│   ├── name/      # Gene name added
│   ├── split/    # Split by chromosome
│   ├── split-uni/  # UniProt added
│   ├── split-uniprot/ # Protein IDs added
│   ├── joined/   # Joined CSVs
├── SAMPLEID-snv.csv  # Final output (example)
