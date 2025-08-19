#!/usr/bin/env python3
"""
validate_build.py

Dry-run validator for the mut-sig-aa-model repo.

What it validates:
  1) Discovers Python scripts under the repo root.
  2) Compiles each with py_compile (syntax check).
  3) Runs --help for key entrypoints (should exit 0 quickly).
  4) Safe import guard for entrypoints (import-time errors).
  5) Optional smoke: run run_sample_pipeline.py with --dry-run.

All outputs are written to --out_dir (default: ./validate_results):
  - validate.log        : human-readable summary
  - report.json         : machine-readable results
  - *_stdout.txt / *_stderr.txt per check (captured subprocess output)

Usage:
  python validate_build.py
  python validate_build.py --out_dir /tmp/val --smoke --print-cmds
  python validate_build.py --python /usr/bin/python3.11 --timeout 60
"""

from __future__ import annotations

import argparse
import json
import os
import shlex
import subprocess
import sys
import time
import traceback
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Dict, List, Optional

# ---- Defaults you can tweak ----
DEFAULT_ENTRYPOINTS = [
    "run_sample_pipeline.py",
    "dna_preprocessor.py",
    "make_manifest.py",
    "multiomic_integration.py",
    "protein_preprocessor.py",  # silently skipped if missing
]

DEFAULT_EXCLUDE_GLOBS = [
    "venv/**", ".venv/**", ".git/**", "__pycache__/**", "build/**", "dist/**",
    "validate_build.py",  # exclude self
]


@dataclass
class CheckResult:
    path: str        # file path or command string
    kind: str        # py_compile | --help | import_guard | smoke_*
    ok: bool
    msg: str = ""
    elapsed: float = 0.0
    returncode: Optional[int] = None
    stdout: Optional[str] = None
    stderr: Optional[str] = None


# ----------------- File helpers -----------------
def ensure_out_dir(d: Path) -> Path:
    d.mkdir(parents=True, exist_ok=True)
    return d

def safe_name_for_files(p: str) -> str:
    # turn "path/to/script.py" + kind into a safe base name
    return p.replace(os.sep, "_").replace(" ", "_").replace(":", "_")

def _print_failures_to_terminal(results: List["CheckResult"], max_stderr_lines: int = 8) -> None:
    """Pretty-print failing checks to the terminal."""
    failures = [r for r in results if not r.ok]
    if not failures:
        return
    print("\nFailures:")
    for r in failures:
        print(f" - {r.kind:18} :: {r.path}")
        if r.msg:
            print(f"    -> {r.msg}")
        if r.stderr:
            lines = r.stderr.strip().splitlines()
            if lines:
                print("    stderr (truncated):")
                for ln in lines[:max_stderr_lines]:
                    print("      " + ln)
                if len(lines) > max_stderr_lines:
                    print("      ...")
        elif r.stdout:
            lines = r.stdout.strip().splitlines()
            if lines:
                print("    stdout (truncated):")
                for ln in lines[:max_stderr_lines]:
                    print("      " + ln)
                if len(lines) > max_stderr_lines:
                    print("      ...")


# ----------------- Subprocess utilities -----------------
def run_cmd(cmd: List[str], timeout: int, env: Dict[str, str], print_cmds: bool) -> subprocess.CompletedProcess:
    if print_cmds:
        print("[cmd]", " ".join(shlex.quote(c) for c in cmd))
    return subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        timeout=timeout,
        env=env,
        text=True,
    )


# ----------------- Individual checks -----------------
def py_compile_check(py: str, filepath: Path, timeout: int, env: Dict[str, str], print_cmds: bool) -> CheckResult:
    t0 = time.time()
    try:
        code = f"import py_compile; py_compile.compile(r'''{filepath}''', doraise=True)"
        proc = run_cmd([py, "-c", code], timeout=timeout, env=env, print_cmds=print_cmds)
        ok = proc.returncode == 0
        msg = "syntax ok" if ok else "py_compile failed"
        return CheckResult(str(filepath), "py_compile", ok, msg, time.time()-t0, proc.returncode, proc.stdout, proc.stderr)
    except Exception as e:
        return CheckResult(str(filepath), "py_compile", False, f"exception: {e}", time.time()-t0, None, None, traceback.format_exc())

def help_check(py: str, script: Path, timeout: int, env: Dict[str, str], print_cmds: bool) -> CheckResult:
    t0 = time.time()
    try:
        proc = run_cmd([py, str(script), "--help"], timeout=timeout, env=env, print_cmds=print_cmds)
        ok = proc.returncode == 0
        msg = "help ok" if ok else "help failed"
        return CheckResult(str(script), "--help", ok, msg, time.time()-t0, proc.returncode, proc.stdout, proc.stderr)
    except Exception as e:
        return CheckResult(str(script), "--help", False, f"exception: {e}", time.time()-t0, None, None, traceback.format_exc())

def import_guard_check(py: str, script: Path, timeout: int, env: Dict[str, str], print_cmds: bool) -> CheckResult:
    """
    Catch import-time errors without running heavy code.
    Run the entrypoint with argv=['script.py','--help'] via runpy.
    """
    t0 = time.time()
    code = (
        "import sys, runpy;"
        f"sys.argv=['{script.name}','--help'];"
        f"runpy.run_path(r'''{script}''', run_name='__main__')"
    )
    try:
        proc = run_cmd([py, "-c", code], timeout=timeout, env=env, print_cmds=print_cmds)
        ok = proc.returncode == 0
        msg = "import+help ok" if ok else "import+help failed"
        return CheckResult(str(script), "import_guard", ok, msg, time.time()-t0, proc.returncode, proc.stdout, proc.stderr)
    except Exception as e:
        return CheckResult(str(script), "import_guard", False, f"exception: {e}", time.time()-t0, None, None, traceback.format_exc())

def shell_check(cmd: List[str], timeout: int, env: Dict[str, str], print_cmds: bool, kind: str, path_label: str) -> CheckResult:
    t0 = time.time()
    try:
        proc = run_cmd(cmd, timeout=timeout, env=env, print_cmds=print_cmds)
        ok = proc.returncode == 0
        msg = "ok" if ok else f"failed (rc={proc.returncode})"
        return CheckResult(path_label, kind, ok, msg, time.time()-t0, proc.returncode, proc.stdout, proc.stderr)
    except Exception as e:
        return CheckResult(path_label, kind, False, f"exception: {e}", time.time()-t0, None, None, traceback.format_exc())


# ----------------- Discovery -----------------
def discover_scripts(root: Path, exclude_globs: List[str]) -> List[Path]:
    all_py = list(root.rglob("*.py"))
    # Expand excludes
    exclude_set = set()
    for pat in exclude_globs:
        for path in root.glob(pat):
            exclude_set.add(path)
    out = [p for p in all_py if p.is_file() and p not in exclude_set]
    return sorted(out)


# ----------------- Validation Orchestration -----------------
def validate_repo(repo_root: Path,
                  python: str,
                  timeout: int,
                  print_cmds: bool,
                  fail_fast: bool,
                  expected_entries: List[str],
                  exclude_globs: List[str],
                  run_smoke: bool,
                  env_overrides: Dict[str, str],
                  out_dir: Path) -> List[CheckResult]:

    results: List[CheckResult] = []
    env = os.environ.copy()
    env.update(env_overrides)

    # 1) Discover & compile
    scripts = discover_scripts(repo_root, exclude_globs)
    for f in scripts:
        res = py_compile_check(python, f, timeout, env, print_cmds)
        results.append(res)
        _capture_outputs(out_dir, res)
        if fail_fast and not res.ok:
            return results

    # 2) entrypoint --help
    for name in expected_entries:
        script = repo_root / name
        if not script.exists():
            continue
        res = help_check(python, script, timeout, env, print_cmds)
        results.append(res)
        _capture_outputs(out_dir, res, label=f"{script.name}_help")
        if fail_fast and not res.ok:
            return results

    # 3) import guard
    for name in expected_entries:
        script = repo_root / name
        if not script.exists():
            continue
        res = import_guard_check(python, script, timeout, env, print_cmds)
        results.append(res)
        _capture_outputs(out_dir, res, label=f"{script.name}_import")
        if fail_fast and not res.ok:
            return results

    # 4) optional smoke of run_sample_pipeline.py --dry-run
    if run_smoke:
        rsp = repo_root / "run_sample_pipeline.py"
        kind = "smoke_run_sample_pipeline"
        if rsp.exists():
            cmd = [python, str(rsp), "--dry-run"]
            res = shell_check(cmd, timeout, env, print_cmds, kind, str(rsp))
        else:
            res = CheckResult(str(rsp), kind, True, "skipped (not present)")
        results.append(res)
        _capture_outputs(out_dir, res, label="run_sample_pipeline_smoke")

    return results


# ----------------- Output capture & reports -----------------
def _capture_outputs(out_dir: Path, res: CheckResult, label: Optional[str] = None) -> None:
    """
    Save stdout/stderr to files under out_dir for each check.
    """
    base = label or f"{safe_name_for_files(res.path)}_{res.kind}"
    if res.stdout:
        (out_dir / f"{base}_stdout.txt").write_text(res.stdout, encoding="utf-8", errors="ignore")
    if res.stderr:
        (out_dir / f"{base}_stderr.txt").write_text(res.stderr, encoding="utf-8", errors="ignore")

def write_log(out_dir: Path, results: List[CheckResult], all_ok: bool = True) -> None:
    """
    Write a human-readable summary to validate.log and append a final banner.
    """
    log_fp = out_dir / "validate.log"
    lines = []
    lines.append("=== validate_build summary ===")
    for r in results:
        status = "OK " if r.ok else "FAIL"
        lines.append(f"[{status}] {r.kind:18} :: {r.path}")
        if r.msg:
            lines.append(f"  -> {r.msg}")
        if r.returncode is not None:
            lines.append(f"  rc={r.returncode}  elapsed={r.elapsed:.2f}s")
    lines.append("")  # blank line before banner

    # Append banner
    if all_ok:
        lines += [
            "==============================",
            "   ✅ BUILD VALIDATED ✅   ",
            "==============================",
        ]
    else:
        lines += [
            "==============================",
            "   ❌ BUILD FAILED ❌   ",
            "==============================",
        ]

    log_fp.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_json(out_dir: Path, results: List[CheckResult], json_name: str = "report.json") -> None:
    payload = {"results": [asdict(r) for r in results]}
    (out_dir / json_name).write_text(json.dumps(payload, indent=2), encoding="utf-8")


# ----------------- CLI -----------------
def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="Dry-run validator for mut-sig-aa-model (repo-wide)")
    ap.add_argument("--repo-root", default=".", help="Repository root (default: .)")
    ap.add_argument("--python", default=sys.executable, help="Python interpreter to use")
    ap.add_argument("--timeout", type=int, default=25, help="Per-subprocess timeout seconds (default: 25)")
    ap.add_argument("--print-cmds", action="store_true", help="Echo subprocess commands")
    ap.add_argument("--fail-fast", action="store_true", help="Stop on first failure")
    ap.add_argument("--smoke", action="store_true", help="Also run run_sample_pipeline.py --dry-run")
    ap.add_argument("--expected-entries", nargs="*", default=DEFAULT_ENTRYPOINTS,
                    help="Entrypoint scripts to check with --help and import guard")
    ap.add_argument("--exclude-glob", dest="exclude_globs", action="append", default=[],
                    help="Glob to exclude (repeatable). Defaults include venv, .git, __pycache__, self")
    ap.add_argument("--env", dest="env_overrides", action="append", default=[],
                    help="Environment override VAR=VAL (repeatable)")
    ap.add_argument("--out_dir", default="validate_results",
                    help="Directory to write logs, report.json, and captured outputs (default: ./validate_results)")
    return ap.parse_args()


def main():
    args = parse_args()
    repo_root = Path(args.repo_root).resolve()
    out_dir = ensure_out_dir(Path(args.out_dir).resolve())

    # build env overrides
    env_overrides: Dict[str, str] = {}
    for item in args.env_overrides:
        if "=" in item:
            k, v = item.split("=", 1)
            env_overrides[k] = v
        else:
            print(f"[warn] ignoring malformed --env {item!r} (expected VAR=VAL)")

    results = validate_repo(
        repo_root=repo_root,
        python=args.python,
        timeout=args.timeout,
        print_cmds=args.print_cmds,
        fail_fast=args.fail_fast,
        expected_entries=args.expected_entries,
        exclude_globs=DEFAULT_EXCLUDE_GLOBS + list(args.exclude_globs or []),
        run_smoke=args.smoke,
        env_overrides=env_overrides,
        out_dir=out_dir,
    )

    # Write artifacts (log + JSON) and final banner
    all_ok = all(r.ok for r in results)
    write_log(out_dir, results, all_ok=all_ok)
    write_json(out_dir, results, json_name="report.json")

    rc = 0 if all_ok else 2
    print(f"\n[validate_build] wrote results to: {out_dir}")
    if all_ok:
        print("\n==============================")
        print("   ✅ BUILD VALIDATED ✅   ")
        print("==============================")
    else:
        print("\n==============================")
        print("   ❌ BUILD FAILED ❌   ")
        print("==============================")
        _print_failures_to_terminal(results)

    print(f"[validate_build] exit code: {rc}")
    sys.exit(rc)


if __name__ == "__main__":
    main()

