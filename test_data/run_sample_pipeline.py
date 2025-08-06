#!/usr/bin/env python3
import subprocess
import os

# Base directories
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
SAMPLE_DATA = os.path.join(BASE_DIR, "test_data")
OBSERVED = os.path.join(BASE_DIR, "observed")
COMPARISON = os.path.join(BASE_DIR, "comparison")
RESULTS_DIR = os.path.join(BASE_DIR, "results")

os.makedirs(RESULTS_DIR, exist_ok=True)

# Paths to repo scripts (assuming this file is in the repo root alongside scripts)
SCRIPTS = {
    "data_prep_summary": "data_prep_summary.py",
    "protein_preprocessor": "protein_preprocessor.py",
    "multiomic_mapper": "multiomic_mapper.py",
    "multiomic_integration": "multiomic_integration.py",
    "multiomic_analysis": "multiomic_analysis.py",
    "signature_modeler": "signature_modeler.py",
    "comparative_analysis": "comparative_analysis.py",
}

# Commands for each step
commands = [
    [
        "python", SCRIPTS["data_prep_summary"],
        os.path.join(SAMPLE_DATA, "mutect"),
        os.path.join(SAMPLE_DATA, "mutect_manifest.txt"),
        os.path.join(RESULTS_DIR, "data_prep"),
        "1",
        "--step", "all"
    ],
    [
        "python", SCRIPTS["protein_preprocessor"],
        "--folder", os.path.join(SAMPLE_DATA, "psm"),
        "--manifest", os.path.join(SAMPLE_DATA, "psm_manifest.txt"),
        "--out", os.path.join(RESULTS_DIR, "psm"),
        "--channel", "Expression",
        "--step", "all"
    ],
    [
        "python", SCRIPTS["multiomic_mapper"],
        "--folder", os.path.join(RESULTS_DIR, "data_prep"),
        "--manifest", os.path.join(SAMPLE_DATA, "mutect_manifest.txt"),
        "--out", os.path.join(RESULTS_DIR, "mapper"),
        "--ref", os.path.join(BASE_DIR, "reference"),
        "--step", "all"
    ],
    [
        "python", SCRIPTS["multiomic_integration"],
        "--manifest", os.path.join(SAMPLE_DATA, "mutect_manifest.txt"),
        "--input_var_dir", os.path.join(RESULTS_DIR, "mapper"),
        "--input_rna_dir", os.path.join(SAMPLE_DATA, "rna"),
        "--rna_manifest", os.path.join(SAMPLE_DATA, "rna_manifest.txt"),
        "--input_ch3_dir", os.path.join(SAMPLE_DATA, "ch3"),
        "--ch3_manifest", os.path.join(SAMPLE_DATA, "ch3_manifest.txt"),
        "--input_protein_dir", os.path.join(RESULTS_DIR, "psm"),
        "--input_cn_dir", os.path.join(SAMPLE_DATA, "cn"),
        "--cn_manifest", os.path.join(SAMPLE_DATA, "cn_manifest.txt"),
        "--out", os.path.join(RESULTS_DIR, "integration")
    ],
    [
        "python", SCRIPTS["multiomic_analysis"],
        "--input_dir", os.path.join(RESULTS_DIR, "integration"),
        "--manifest", os.path.join(SAMPLE_DATA, "mutect_manifest.txt"),
        "--out_dir", os.path.join(RESULTS_DIR, "analysis"),
        "--step", "all"
    ],
    [
        "python", SCRIPTS["signature_modeler"],
        "--observed_dir", os.path.join(RESULTS_DIR, "data_prep"),
        "--manifest", os.path.join(SAMPLE_DATA, "mutect_manifest.txt"),
        "--out", os.path.join(RESULTS_DIR, "comparison"),
        "--step", "all"
    ],
    [
        "python", SCRIPTS["comparative_analysis"],
        "--observed_dir", OBSERVED,
        "--comparison_dir", COMPARISON,
        "--manifest", os.path.join(SAMPLE_DATA, "mutect_manifest.txt"),
        "--out_dir", os.path.join(RESULTS_DIR, "comparative"),
        "--step", "all"
    ]
]

# Run steps with error handling
for i, cmd in enumerate(commands, start=1):
    step_name = os.path.basename(cmd[1])
    print(f"\n=== Step {i}: {step_name} ===")
    try:
        subprocess.run(cmd, check=True)
        print(f"Step {i} completed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Step {i} FAILED with error: {e.stderr if e.stderr else e}")

print("\nPipeline test run complete. Results in:", RESULTS_DIR)
