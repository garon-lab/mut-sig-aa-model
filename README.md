# MUT-SIG-AA-MODEL
Uses publicly available genomic data to create multiomic structures of amino acid substitutions from non-synonymous mutations with ability to compare amino acid substitution matrices and model expected output from mutational signature. 

You will need downloaded VEP-annotated mutect, RNA (RNASeq), CH3 (methylation), PSM (protein), and copy number (CN) files with their respective manifests in standard GDC format (https://proteomic.datacommons.cancer.gov/pdc/). Standard references are available for dowload with this script and otherwise should include esng gene symbols, gene symbol names, uni-ensg, uni-np, and methylation probe IDs matched to chromosome (CHR) and gene name.

Python 3.8+ is required.

Recommended order of scripts:
1. DNA (VEP) preprocessor
2. Protein (PSM) preprocesser
3. Make manifest
4. Multiomic data integration
5. Multiomic analysis
6. Signature modeler (optional)
7. Comparative analysis (optional)

Additional information regarding each script is below and can be accessed with the --help flag.


# QUICK START
To check install, the run_sample_pipeline script will download dependencies and use test_data may be used to ensure all components are functioning as expected.

Usage example:
python run_sample_pipeline.py

# DNA PREPROCESSOR

This script pre-processes VEP files to select for any variants of interest. The example in our manuscript is selecting for single nucleotide variants (SNVs) and sningle nucleotide polymorphisms (SNPs), but can be adapted as needed by the user.



# PROTEIN PREPROCESSOR

This script processes protein PSM files by:
1. Splitting a sample manifest into channels.
2. Creating index folders for output.
3. Filtering and reorganizing raw data by TMT channel.
4. Rejoining processed files into complete datasets named by sample-id.

Dependencies: pandas, numpy, argparse, pathlib

Usage:
    python protein_preprocessor.py \
        --folder <input directory> \
        --manifest <manifest file> \
        --out_dir <output directory> \
        --channel <TMT channel (e.g., 126,127N)> \
        [--step all|channels|index|prep|join]

Arguments:
    --folder    Directory with raw input files
    --manifest  Tab-delimited manifest file listing sample metadata
    --out_dir   Output directory to store results
    --channel   TMT channel identifier (e.g., 126, 127N, etc.)
    --step      Pipeline step to run: all (default), channels, index, prep, or join


# MAKE MANIFEST

This script organizes metafiles (manifests) that will be used to integrate multiomic information for calling during multiomic integration. Outputs are tab-separated TSVs.

Dependencies: pandas

Project layout expected (for --build main)

PROJECT_ROOT/
   dna/        (optional)         case-id.csv                       #VEP-annotated mutect files
               gdc-dna.tsv        + dna/{File ID}/{File Name}
   rna/        gdc-rna.tsv        + rna/{File ID}/{File Name}       #RNA-seq count data
   ch3/        gdc-ch3.tsv        + ch3/{File ID}/{File Name}       #Methylation data in GDC format
   cn/         gdc-cn.tsv         + cn/{File ID}/{File Name}        #Copy Number data in GDC format
   protein/    (optional)         case-id.csv                       #PSM reads (NP, SEQ, EV, INT)

Usage (recommended):
   python make_manifest.py \
        --build-main \
        --project_root <project directory> \
        --main_out <output directory>

Usage (single-type):
   python make_manifest.py \
        --gdc-manifest <manifest file> \
        --out_dir <output directory> \
        --type <rna|ch3|cn> \
        --folder-id-col <# (e.g., 0)> \
        --file-name-col <# (e.g., 1)> \
        --case-id-col <# (e.g, 5)> 

Arguments:
(Recommended Use)
   --build-main        Flag that creates unified manifest with paths constructed like {PROJECT_ROOT}/{modality}/{File ID}/{File Name}
   --project_root      Directory that has all project-related files (dna, rna, ch3, cn, protein)
   --main_out          Directory that will contain all output files
   
(Single Use)   
   --gdc-manifest      Flag to use if opting for creation of single manifest as opposed to build-main
   --out_dir           Directory that will contain output file
   --type              Used for naming output manifest (e.g., rna_manifest.tsv)
   --folder-id-column  Column number that identifies path (downloaded code path; e.g, 123abc)
   --file-name-col     Column number that identifies file name (downloaded code + file name; e.g., 987zyx.wxs.MuTect2.somatic_annotation.vcf.gz)
   --case-id-col       Column number that identifies case-id (e.g., C3N-001)

   
# MULTIOMIC INTEGRATION

This pipeline integrates multiple omics layers (DNA, RNA, methylation, protein, and copy number)
for each case-id listed in a manifest. Can be used with an existing unifed manifest built with make_manifest.py or separate manifestsOnly the final protein files (with copy number added) are retained
in the output directory. Note protein files must be preprocessed and in the format {Case-ID}.csv (see protein_preprocessor.py)

Dependencies: csv, argparse, shutil, pandas, numpy

Recommended Usage:


Optional Usage: 
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


# MULTIOMIC ANALYSIS

This script provides downstream analysis for integrated multiomic data.

Features:
1. Statistical summaries of each omic layer across samples
2. Plot generation (e.g., heatmaps, scatter plots)
3. Clustering and correlation analysis
4. Protein-expression centered summarization (protein_only, no_protein, single_entry per gene)

Dependencies: argparse, logging, pathlib, pandas, numpy, matplotlib, scipy

Usage:
    python multiomic_analysis.py \
        --input_dir <directory of integrated CSVs> \
        --manifest <manifest file> \
        --out_dir <output directory> \
        [--step stats|plots|cluster|protein_only|no_protein|single_entry|all]

Arguments:
    --input_dir  Directory containing per-sample integrated CSV files
    --manifest   Tab-delimited manifest file listing sample IDs in the first column
    --out_dir    Directory to write analysis outputs
    --step       Analysis step to run: stats, plots, cluster, protein_only (selects for genes with protein expression), no_protein (selects for genes without protein expression), single_entry (compresses each gene information into a single row selecting for highest E value and averaging beta values), or all (default: all)


# SIGNATURE MODELER

This repository provides a standalone script to model expected amino acid substitution profiles from mutational signature context proportions.

Features:
1. Provides expected amino acid substitution matrices (21x21) using embedded base matrices.
2. Performs row-wise normalization to predefined amino-acid frequency targets based on GRCh38 amino-acid frequency targets.
3. Flattens vectors for downstream analysis (expected_vectors.csv).
4. Optional heatmap visualization of expected substitution profiles (expected_heatmap.png).

Dependencies: pandas, numpy, matplotlib, seaborn

Usage:
   python signature_modeling.py \
    --signature_vector contexts.csv \
    --out_dir results/ \
    [--step model|heatmap|both] \
    [--log_level DEBUG|INFO|WARNING|ERROR]

Arguments:
    --signature_vector  Path to CSV file with columns [ID, AC, AG, AT, CA, CG, CT, GA, GC, GT, TA, TC, TG]
    --out_dir           Directory where outputs will be written (created if necessary)
    --step              Pipeline step(s) to run: model, heatmap, both (default both)
    --log_level         Logging verbosity (default INFO)


# COMPARATIVE ANALYSIS

This script summarizes observed amino acid variant counts per sample, compares amino acid substitution profiles from observed or modeled data (e.g., from signature modeler), and visualizes results as heatmaps.

Features:
1. Aggregation of observed amino acid variants into a single vector per sample.
2. Calculates cosine similarity between observed variant vectors and/or modeled vectors.
3. Generates heatmaps of similarity matrices or per-sample counts.

Dependencies: pip, install, pandas, numpy, matplotlib, seaborn

Usage:
  python comparative_analysis.py \
    --observed_dir <directory of observed AA matrix CSVs> \
    --comparison_dir <directory of comparison AA matrix CSVs> \
    --manifest <manifest file> \
    --out_dir <output directory> \
    [--vector_file <observed summary CSV>] \
    [--step summarize|compare|heatmap|single-file|all]

Arguments:
    --observed_dir      Directory of observed amino acid csvs, names {sample-id}.csv
    --comparison_dir    Directory of comparison amino acid csvs,  named {sample-id}.csv
    --mainfest          Tab-deliminted manifest listing sample IDs (first column)
    --out_dir           Output directory for csvs and plots
    --vector_file       Path to save/load summary csvs (default: out_dir/observed_summary.csv)
    --step              Pipeline step(s) to run: summarize, compare, heatmap, single-file, or all (default: all)

Output structure
results/
├── observed_summary.csv        # Summarized observed counts
├── similarity_matrix.csv       # Cosine similarity scores
├── heatmap.png                  # Similarity matrix heatmap
├── aa-count.png                 # Per-sample count heatmap
└── aa-proportion.png            # Per-sample proportion heatmap
