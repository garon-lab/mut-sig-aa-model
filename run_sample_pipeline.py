#!/usr/bin/env python3
import os
import sys
import subprocess
import time
import urllib.request
import zipfile
from pathlib import Path
import importlib.util

# -----------------------
# Config
# -----------------------
ZIP_URL = "https://github.com/garon-lab/mut-sig-aa-model/raw/main/test_data.zip"
ZIP_NAME = "test_data.zip"
DATA_DIR = "test_data"
RESULTS_DIR = "results"
SUMMARY_FILE = os.path.join(RESULTS_DIR, "summary.txt")
REPO_ROOT = Path(__file__).parent.resolve()

REQUIRED_SCRIPTS = [
    "data_prep_summary.py",
    "protein_preprocessor.py",
    "multiomic_mapper.py",
    "multiomic_integration.py",
    "multiomic_analysis.py",
    "signature_modeler.py",
    "comparative_analysis.py"
]

REQUIRED_PKGS = ["pandas", "numpy", "matplotlib", "seaborn"]

# ANSI Colors
GREEN = "\033[92m"
RED = "\033[91m"
RESET = "\033[0m"

# -----------------------
# Helpers
# -----------------------
def install_dependencies():
    for pkg in REQUIRED_PKGS:
        if importlib.util.find_spec(pkg) is None:
            print(f"Installing missing dependency: {pkg}")
            subprocess.check_call([sys.executable, "-m", "pip", "install", pkg])
        else:
            print(f"Dependency OK: {pkg}")

def download_data():
    if not Path(ZIP_NAME).exists():
        print(f"Downloading {ZIP_NAME} from {ZIP_URL} ...")
        urllib.request.urlretrieve(ZIP_URL, ZIP_NAME)
        print("Download complete.")

def extract_data():
    if not Path(DATA_DIR).exists():
        print(f"Extracting {ZIP_NAME} ...")
        with zipfile.ZipFile(ZIP_NAME, "r") as zip_ref:
            zip_ref.extractall(DATA_DIR)
        print("Extraction complete.")

def ensure_gitignore():
    gitignore_path = REPO_ROOT / ".gitignore"
    if gitignore_path.exists():
        content = gitignore_path.read_text().splitlines()
    else:
        content = []
    if f"{DATA_DIR}/" not in content:
        content.append(f"{DATA_DIR}/")
        gitignore_path.write_text("\n".join(content))
        print(f"Added {DATA_DIR}/ to .gitignore")

def run_step(script, args):
    cmd = [sys.executable, str(REPO_ROOT / script)] + args
    start = time.time()
    try:
        subprocess.run(cmd, check=True)
        duration = time.time() - start
        return ("PASS", duration, "Completed successfully")
    except subprocess.CalledProcessError as e:
        duration = time.time() - start
        return ("FAIL", duration, e.stderr.decode() if e.stderr else str(e))

def main():
    os.makedirs(RESULTS_DIR, exist_ok=True)
    install_dependencies()
    download_data()
    extract_data()
    ensure_gitignore()

    # Detect scripts
    scripts_found = [s for s in REQUIRED_SCRIPTS if (REPO_ROOT / s).exists()]
    if len(scripts_found) < len(REQUIRED_SCRIPTS):
        print(f"Warning: Missing scripts: {set(REQUIRED_SCRIPTS) - set(scripts_found)}")

    # Step arguments
    step_args = {
        "data_prep_summary.py": [
            f"{DATA_DIR}/test_data/test_mutect",
            f"{DATA_DIR}/test_data/mutect_manifest.tsv",
            f"{RESULTS_DIR}/data_prep",
            "2", "--step", "all"
        ],
        "protein_preprocessor.py": [
            "--folder", f"{DATA_DIR}/test_data/psm",
            "--manifest", f"{DATA_DIR}/test_data/psm_manifest.txt",
            "--out", f"{RESULTS_DIR}/psm",
            "--channel", "Expression",
            "--step", "all"
        ],
        "multiomic_mapper.py": [
            "--folder", f"{RESULTS_DIR}/data_prep",
            "--manifest", f"{DATA_DIR}/test_data/mutect_manifest.txt",
            "--out", f"{RESULTS_DIR}/mapper",
            "--ref", f"{REPO_ROOT}/reference",
            "--step", "all"
        ],
        "multiomic_integration.py": [
            "--manifest", f"{DATA_DIR}/test_data/mutect_manifest.txt",
            "--input_var_dir", f"{RESULTS_DIR}/mapper",
            "--input_rna_dir", f"{DATA_DIR}/test_data/rna",
            "--rna_manifest", f"{DATA_DIR}/test_data/rna_manifest.txt",
            "--input_ch3_dir", f"{DATA_DIR}/test_data/ch3",
            "--ch3_manifest", f"{DATA_DIR}/test_data/ch3_manifest.txt",
            "--input_protein_dir", f"{RESULTS_DIR}/psm",
            "--input_cn_dir", f"{DATA_DIR}/test_data/cn",
            "--cn_manifest", f"{DATA_DIR}/test_data/cn_manifest.txt",
            "--out", f"{RESULTS_DIR}/integration"
        ],
        "multiomic_analysis.py": [
            "--input_dir", f"{RESULTS_DIR}/integration",
            "--manifest", f"{DATA_DIR}/test_data/mutect_manifest.txt",
            "--out_dir", f"{RESULTS_DIR}/analysis",
            "--step", "all"
        ],
        "signature_modeler.py": [
            "--observed_dir", f"{RESULTS_DIR}/data_prep",
            "--manifest", f"{DATA_DIR}/test_data/mutect_manifest.txt",
            "--out", f"{RESULTS_DIR}/comparison",
            "--step", "all"
        ],
        "comparative_analysis.py": [
            "--observed_dir", f"{DATA_DIR}/observed",
            "--comparison_dir", f"{DATA_DIR}/comparison",
            "--manifest", f"{DATA_DIR}/test_data/mutect_manifest.txt",
            "--out_dir", f"{RESULTS_DIR}/comparative",
            "--step", "all"
        ]
    }

    summary = []
    for script in scripts_found:
        print(f"\n=== Running {script} ===")
        status, duration, message = run_step(script, step_args[script])
        color_status = f"{GREEN}{status}{RESET}" if status == "PASS" else f"{RED}{status}{RESET}"
        print(f"{script} - {color_status} ({duration:.2f}s) - {message}")
        summary.append((script, status, duration, message))

    # Print and save summary table
    header = f"{'Step':<30} {'Status':<8} {'Time (s)':<10} Message"
    sep = "-" * len(header)
    print("\n" + header)
    print(sep)
    with open(SUMMARY_FILE, "w") as f:
        f.write(header + "\n" + sep + "\n")
        for step, status, duration, message in summary:
            status_col = status
            color_status = f"{GREEN}{status}{RESET}" if status == "PASS" else f"{RED}{status}{RESET}"
            line = f"{step:<30} {color_status:<8} {duration:<10.2f} {message}"
            print(line)
            f.write(f"{step:<30} {status_col:<8} {duration:<10.2f} {message}\n")

    print(f"\nSummary saved to {SUMMARY_FILE}")

if __name__ == "__main__":
    main()
