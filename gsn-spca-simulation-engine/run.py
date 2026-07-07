from __future__ import annotations

import argparse
import configparser
import json
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parent
DEFAULT_JSON = ROOT / "configs" / "default_config.json"
DEFAULT_SETTINGS = ROOT / "configs" / "run_settings.txt"
RUNTIME_CONFIG_DIR = ROOT / "outputs" / "runtime_configs"


def parse_bool(value: str) -> bool:
    return value.strip().lower() in {"1", "true", "yes", "y", "on"}


def parse_optional_int(value: str) -> int | None:
    value = value.strip()
    if not value:
        return None
    return int(value)


def parse_param_value(key: str, raw: str) -> Any:
    raw = raw.strip()
    if raw.lower() == "auto_from_summary":
        return "auto_from_summary"
    if key in {"k", "niter", "num_init"}:
        return int(raw)
    return float(raw)


def apply_settings(base_config: dict[str, Any], settings_path: Path) -> tuple[dict[str, Any], dict[str, str]]:
    parser = configparser.ConfigParser()
    parser.optionxform = str
    parser.read(settings_path, encoding="utf-8-sig")

    run_settings = {
        "mode": parser.get("run", "mode", fallback="collect").strip(),
        "methods": parser.get("run", "methods", fallback="all").strip(),
        "parameter_sets": parser.get("run", "parameter_sets", fallback="default").strip(),
        "run_id": parser.get("run", "run_id", fallback="").strip(),
        "prepared_inputs": parser.get("run", "prepared_inputs", fallback="").strip(),
    }
    base_config["selected_parameter_sets"] = run_settings["parameter_sets"]

    if parser.has_section("paths"):
        path_keys = [
            "source_root",
            "rscript",
            "reticulate_python",
            "parameter_set_input_root",
            "parameter_set_collect_summary",
        ]
        for key in path_keys:
            value = parser.get("paths", key, fallback="").strip()
            if value:
                base_config[key] = value

    base_config.setdefault("parameter_sets_csv", str(ROOT / "configs" / "parameter_sets.csv"))

    if parser.has_section("seeds"):
        seeds = base_config.setdefault("seeds", {})
        for key, value in parser.items("seeds"):
            value = value.strip()
            if value:
                try:
                    seeds[key] = int(value)
                except ValueError:
                    seeds[key] = value

    for method in base_config.get("methods", {}):
        section = f"method.{method}"
        if not parser.has_section(section):
            continue
        method_cfg = base_config["methods"][method]
        params = method_cfg.setdefault("params", {})
        for key in ["k", "k_group", "we", "t", "niter", "err", "num_init"]:
            if parser.has_option(section, key):
                raw = parser.get(section, key).strip()
                if raw:
                    params[key] = parse_param_value(key, raw)
        if parser.has_option(section, "edge_mode"):
            edge_mode = parser.get(section, "edge_mode").strip()
            if edge_mode:
                method_cfg["edge_mode"] = edge_mode
        if parser.has_option(section, "hardcoded_edge_n"):
            method_cfg["hardcoded_edge_n"] = parse_optional_int(parser.get(section, "hardcoded_edge_n"))
        if method == "GSN-SPCA":
            espca_lee = method_cfg.setdefault("espca_lee", {})
            if parser.has_option(section, "espca_lee_enabled"):
                espca_lee["enabled"] = parse_bool(parser.get(section, "espca_lee_enabled"))
            if parser.has_option(section, "espca_lee_N"):
                raw = parser.get(section, "espca_lee_N").strip()
                if raw:
                    espca_lee["N"] = float(raw)
            if parser.has_option(section, "espca_lee_Q"):
                raw = parser.get(section, "espca_lee_Q").strip()
                if raw:
                    espca_lee["Q"] = float(raw)

    return base_config, run_settings


def build_runtime_config(json_config: Path, settings_path: Path) -> tuple[Path, dict[str, str]]:
    base_config = json.loads(json_config.read_text(encoding="utf-8-sig"))
    config, run_settings = apply_settings(base_config, settings_path)
    RUNTIME_CONFIG_DIR.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    runtime_path = RUNTIME_CONFIG_DIR / f"runtime_config_{stamp}.json"
    runtime_path.write_text(json.dumps(config, ensure_ascii=False, indent=2), encoding="utf-8")
    return runtime_path, run_settings


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Python-first entry for the GSN-SPCA simulation experiment package."
    )
    parser.add_argument("--settings", default=str(DEFAULT_SETTINGS), help="User-editable text settings file.")
    parser.add_argument("--base-config", default=str(DEFAULT_JSON), help="Internal base JSON config.")
    parser.add_argument("--mode", choices=["collect", "run"], default="", help="Override [run] mode.")
    parser.add_argument("--methods", default="", help="Override [run] methods.")
    parser.add_argument("--parameter-sets", default="", help="Override [run] parameter_sets.")
    parser.add_argument("--run-id", default="", help="Optional output run folder name override.")
    parser.add_argument("--prepared-inputs", default="", help="Optional prepared input folder override.")
    parser.add_argument("--check-env", action="store_true", help="Check Python/R/path requirements and exit.")
    args = parser.parse_args()

    if args.check_env:
        check_script = ROOT / "scripts" / "check_environment.py"
        return subprocess.call([sys.executable, str(check_script), "--settings", str(Path(args.settings).resolve())], cwd=ROOT)

    settings_path = Path(args.settings).resolve()
    runtime_config, run_settings = build_runtime_config(Path(args.base_config).resolve(), settings_path)

    mode = args.mode or run_settings["mode"] or "collect"
    methods = args.methods or run_settings["methods"] or "all"
    parameter_sets = args.parameter_sets or run_settings["parameter_sets"] or "default"
    run_id = args.run_id or run_settings["run_id"]
    prepared_inputs = args.prepared_inputs or run_settings["prepared_inputs"]

    cmd = [
        sys.executable,
        str(ROOT / "scripts" / "engine.py"),
        "--config",
        str(runtime_config),
        "--mode",
        mode,
        "--methods",
        methods,
        "--parameter-sets",
        parameter_sets,
    ]
    if run_id:
        cmd.extend(["--run-id", run_id])
    if prepared_inputs:
        cmd.extend(["--prepared-inputs", prepared_inputs])

    print(f"Settings: {settings_path}")
    print(f"Runtime config: {runtime_config}")
    print(f"Mode: {mode}")
    print(f"Parameter sets: {parameter_sets}")
    print(f"Methods: {methods}")
    return subprocess.call(cmd, cwd=ROOT)


if __name__ == "__main__":
    raise SystemExit(main())
