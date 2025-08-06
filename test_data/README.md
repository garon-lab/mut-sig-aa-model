# MUT-SIG-AA-MODEL Sample Data

This folder contains **synthetic example data** for testing the full MUT-SIG-AA-MODEL pipeline.
It is designed for functional testing only and does **not** represent real biological data.

## Directory Structure

```
sample_data/
├── test_data/         # Synthetic CPTAC-style input data
│   ├── mutect/        # VEP-annotated variant files
│   ├── psm/           # Protein expression files
│   ├── rna/           # RNASeq expression files
│   ├── ch3/           # Methylation data files
│   ├── cn/            # Copy number variation files
│   ├── *_manifest.txt # Manifest files listing sample IDs
│
├── observed/          # Dummy observed summary CSVs for comparative_analysis.py
│   ├── SAMPLE1.csv
│   ├── SAMPLE2.csv
│
├── comparison/        # Dummy comparison output CSVs for comparative_analysis.py
│   ├── SAMPLE1.csv
│   ├── SAMPLE2.csv
```

## Usage

You can use these files to run the full pipeline without needing real CPTAC data:

```bash
# Example: Comparative Analysis
python comparative_analysis.py     --observed_dir sample_data/observed     --comparison_dir sample_data/comparison     --manifest sample_data/test_data/mutect_manifest.txt     --out_dir results     --step all
```

## Notes
- File contents are **synthetic** and randomly generated.
- Intended only for testing **script functionality** and **argument parsing**.
- Replace with real data for scientific analyses.
