# MUT-SIG-AA-MODEL

[![Python Version](https://img.shields.io/badge/python-3.8%2B-blue.svg)](https://www.python.org/downloads/release/python-380/)
[![License](https://img.shields.io/badge/license-UCLA-green.svg)](./LICENSE)
[![Dependencies](https://img.shields.io/badge/dependencies-pandas%20%7C%20numpy%20%7C%20matplotlib%20%7C%20seaborn-orange.svg)]()

This repository uses publicly available genomic data to create **multiomic structures of amino acid substitutions** from non-synonymous mutations. It allows comparison of amino acid substitution matrices and modeling of expected output from mutational signatures.  

You will need downloaded **VEP-annotated Mutect**, **RNA (RNA-seq)**, **CH3 (methylation)**, **PSM (protein)**, and **copy number (CN)** files with their respective manifests in **standard GDC format** ([GDC Data Portal](https://proteomic.datacommons.cancer.gov/pdc/)).  

Standard references are available for download with this repository. Otherwise, references should include:  
- **ESNG gene symbols**  
- **Gene symbol names**  
- **Uni-ENSG IDs**  
- **Uni-NP IDs**  
- **Methylation probe IDs matched to chromosome (CHR) and gene name**  

**Requirements:** Python 3.8+  

---

## Table of Contents
- [Recommended Script Order](#recommended-script-order)  
- [Quick Start](#quick-start)  
- [DNA Preprocessor](#dna-preprocessor)  
- [Protein Preprocessor](#protein-preprocessor)  
- [Make Manifest](#make-manifest)  
- [Multiomic Integration](#multiomic-integration)
- [DNA Analysis](#multiomic-analysis)  
- [Multiomic Analysis](#multiomic-analysis)  
- [Signature Modeler](#signature-modeler)  
- [Comparative Analysis](#comparative-analysis)  
- [License](#license)  

---

## Recommended Script Order
1. **DNA Preprocessor** (VEP-annotated Mutect → per-case DNA CSVs)  
2. **Protein Preprocessor** (PSM reads → per-case protein CSVs)  
3. **Make Manifest** (unified or per-omic manifests)  
4. **Multiomic Data Integration**
5. **DNA Analysis**  
6. **Multiomic Analysis**  
7. **Signature Modeler** *(optional)*  
8. **Comparative Analysis** *(optional)*  

Each script provides detailed help with the `--help` flag.  

---

## Quick Start

To verify installation, you can:  

- Run `validate_build.py` to check syntax, entrypoints, and optional smoke tests  
- Run `run_sample_pipeline.py` with the provided test.zip to confirm all components  

```bash
# Lightweight build validation
python validate_build.py --out_dir validate_results

# Optional smoke test
python validate_build.py --smoke

# Full functional test
python run_sample_pipeline.py -out_dir OUT_DIR
```

---

## DNA Preprocessor

Creates per-case DNA CSVs under `dna/{Case-ID}.csv` for use in integration.  

**Features:**  
- Parses VEP-annotated Mutect files  
- Filters for SNVs/SNPs  
- Optional summary, mutational signatures, amino acid substitutions, and matrices  

**Dependencies:** `pandas`

### Usage
```bash
python dna_preprocessor.py   --folder <input directory>   --manifest <gdc-manifest tsv>   --out_dir <output directory>
```

**Optional:**  
```bash
  --make-simplified   --preprocess-mutect   --vcf-folder <vcf directory>   --simplified <case_ids.txt>   --write-signatures --signature-label <label>   --extract-mutations <snp|snv> --write-matrices
```

---

## Protein Preprocessor

Creates per-case protein CSVs under `protein/{Case-ID}.csv`.  

**Features:**  
- Splits manifest into channels  
- Prepares per-channel outputs  
- Joins PSM parts back into per-case datasets  

**Dependencies:** `pandas`, `numpy`, `argparse`, `pathlib`

### Usage
```bash
python protein_preprocessor.py   --folder <input directory>   --manifest <gdc-manifest>   --out_dir <output directory>   --channel <channels>   --step <all|split|prep|join>   --jobs 8
```

---

## Make Manifest

Builds project-wide or per-modality manifests.  

**Dependencies:** `pandas`

### Usage
```bash
python make_manifest.py   --build-main   --project_root <project root>   --out_dir <output directory>
```

*Single-modality example:*  
```bash
python make_manifest.py   --gdc-manifest <manifest file>   --out_dir <output directory>   --type <rna|ch3|cn>   --folder-id-col <#>   --file-name-col <#>   --case-id-col <#>
```

---

## Multiomic Integration

Integrates DNA, RNA, methylation (CH3), protein, and CNV.  
Only final **protein-centered integrated files** are retained.  

**Dependencies:** `csv`, `argparse`, `shutil`, `pandas`, `numpy`

### Usage
```bash
python multiomic_integration.py   --folder <project root>   --manifest <manifest file>   --out_dir <output directory>   
```

**Optional:**  
```bash
  --input_dna_dir DNA_DIR   --input_rna_dir RNA_DIR --rna_manifest RNA_MANIFEST   --input_ch3_dir CH3_DIR --ch3_manifest CH3_MANIFEST   --input_protein_dir PROTEIN_DIR   --input_cn_dir CNV_DIR --cn_manifest CNV_MANIFEST   --skip_rna --skip_ch3 --skip_protein --skip_cn --step <all|dna|rna|ch3|cnv|protein>
```
---

## DNA Analysis

Provides downstream analysis for DNA data.  

**Features:**  
- Case-level summary of SNV and SNP counts  
- Signature composition (12-channel REF->ALT) per case  
- Amino acid substitution extract (from INFO field) 
- Amino acid substitution matrices (21x21) separated by label (SNV vs SNP, or custom labels)

**Dependencies:** `argparse`, `logging`, `pathlib`, `pandas`, `concurrent.futures`

### Usage
```bash
python dna_analysis.py   --in_dir <directory with dna/{Case-ID}.csv>   --simplified <case_ids.txt>   --out_dir <output directory> \
  [--summarize-variants   --write-signatures   --extract-mutations   --write-matrices   --labels Label1 Label2   --max-records N --jobs N]
```
---

## Multiomic Analysis

Provides downstream analysis for integrated data.  

**Features:**  
- Statistical summaries  
- Plot generation (heatmaps, scatter plots)  
- Clustering and correlations  
- Protein-centric summarization  

**Dependencies:** `argparse`, `logging`, `pathlib`, `pandas`, `numpy`, `matplotlib`, `scipy`

### Usage
```bash
python multiomic_analysis.py   --input_dir <integrated CSVs>   --manifest <manifest file>   --out_dir <output directory>   [--step stats|plots|cluster|protein_only|no_protein|single_entry|all]
```

---

## Signature Modeler

Models expected amino acid substitution profiles from mutational signatures.  

**Features:**  
- Generates expected substitution matrices (21×21)  
- Normalizes to GRCh38 amino acid frequencies  
- Outputs vectors and optional heatmaps  

**Dependencies:** `pandas`, `numpy`, `matplotlib`, `seaborn`

### Usage
```bash
python signature_modeling.py   --signature_vector contexts.csv   --out_dir results/   [--step model|heatmap|both]   [--log_level DEBUG|INFO|WARNING|ERROR]
```

---

## Comparative Analysis

Compares observed vs. modeled amino acid substitution profiles.  

**Features:**  
- Aggregates observed variants  
- Computes cosine similarity  
- Generates heatmaps  

**Dependencies:** `pandas`, `numpy`, `matplotlib`, `seaborn`

### Usage
```bash
python comparative_analysis.py   --observed_dir <observed AA CSVs>   --comparison_dir <comparison AA CSVs>   --manifest <manifest file>   --out_dir <output directory>   [--vector_file observed_summary.csv]   [--step summarize|compare|heatmap|single-file|all]
```

---

## License

© UCLA, Edward Garon & Amy Cummings, 2025
