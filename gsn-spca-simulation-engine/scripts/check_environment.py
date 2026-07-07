from __future__ import annotations

import argparse
import configparser
import importlib.util
import json
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SETTINGS = ROOT / "configs" / "run_settings.txt"
DEFAULT_CONFIG = ROOT / "configs" / "default_config.json"


def read_settings(path: Path) -> configparser.ConfigParser:
    parser = configparser.ConfigParser()
    parser.optionxform = str
    parser.read(path, encoding="utf-8-sig")
    return parser


def setting(parser: configparser.ConfigParser, section: str, key: str, fallback: str = "") -> str:
    return parser.get(section, key, fallback=fallback).strip()


def exists(path_text: str) -> bool:
    return bool(path_text) and Path(path_text).exists()


def check_python_module(name: str) -> tuple[str, bool, str]:
    spec = importlib.util.find_spec(name)
    return name, spec is not None, "" if spec else "not installed"


def run_command(cmd: list[str], timeout: int = 30) -> tuple[int, str, str]:
    try:
        proc = subprocess.run(cmd, text=True, capture_output=True, timeout=timeout)
        return proc.returncode, proc.stdout.strip(), proc.stderr.strip()
    except Exception as exc:
        return 999, "", str(exc)


def check_r_package(rscript: str, package: str) -> tuple[str, bool, str]:
    code = f'quit(status = ifelse(requireNamespace("{package}", quietly = TRUE), 0, 1))'
    rc, stdout, stderr = run_command([rscript, "-e", code])
    return package, rc == 0, stderr or stdout


def print_status(ok: bool, label: str, detail: str = "") -> None:
    prefix = "OK" if ok else "MISSING"
    if detail:
        print(f"[{prefix}] {label}: {detail}")
    else:
        print(f"[{prefix}] {label}")


def main() -> int:
    parser = argparse.ArgumentParser(description="Check local environment for the simulation engine package.")
    parser.add_argument("--settings", default=str(DEFAULT_SETTINGS))
    args = parser.parse_args()

    settings_path = Path(args.settings).resolve()
    if not settings_path.exists():
        print_status(False, "settings file", str(settings_path))
        return 2

    settings = read_settings(settings_path)
    run_mode = setting(settings, "run", "mode", "collect")
    rscript = setting(settings, "paths", "rscript")
    reticulate_python = setting(settings, "paths", "reticulate_python")
    source_root = setting(settings, "paths", "source_root")
    parameter_set_input_root = setting(settings, "paths", "parameter_set_input_root")
    parameter_set_collect_summary = setting(settings, "paths", "parameter_set_collect_summary")

    print(f"Settings: {settings_path}")
    print(f"Configured mode: {run_mode}")
    print(f"Current Python: {sys.executable}")
    print(f"Python version: {sys.version.split()[0]}")

    ok_all = True
    py_ok = sys.version_info >= (3, 9)
    print_status(py_ok, "Python >= 3.9")
    ok_all = ok_all and py_ok

    for label, value in [
        ("source_root", source_root),
        ("parameter_set_input_root", parameter_set_input_root),
        ("parameter_set_collect_summary", parameter_set_collect_summary),
    ]:
        ok = exists(value)
        print_status(ok, label, value)
        ok_all = ok_all and ok

    rscript_ok = exists(rscript)
    print_status(rscript_ok, "Rscript path", rscript)
    if run_mode == "run":
        ok_all = ok_all and rscript_ok

    reticulate_ok = exists(reticulate_python)
    print_status(reticulate_ok, "reticulate_python path", reticulate_python)
    if run_mode == "run":
        ok_all = ok_all and reticulate_ok

    print("Python module checks:")
    for module in ["pandas", "sklearn"]:
        name, ok, detail = check_python_module(module)
        print_status(ok, name, detail)
        if run_mode == "run":
            ok_all = ok_all and ok

    if rscript_ok:
        rc, stdout, stderr = run_command([rscript, "--version"])
        print_status(rc == 0, "Rscript executable", stdout or stderr)
        print("R package checks:")
        for package in ["reticulate", "readr", "dplyr"]:
            name, ok, detail = check_r_package(rscript, package)
            print_status(ok, name, detail)
            if run_mode == "run" and package == "reticulate":
                ok_all = ok_all and ok
    else:
        print("R package checks skipped because Rscript was not found.")

    if ok_all:
        print("Environment check passed for the configured mode.")
        return 0
    print("Environment check found missing items. Edit configs/run_settings.txt or install the missing dependencies.")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
