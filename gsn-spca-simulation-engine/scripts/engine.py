from __future__ import annotations

import argparse
import csv
import html
import json
import os
import shutil
import subprocess
from datetime import datetime
from pathlib import Path
from typing import Any


PACKAGE_ROOT = Path(__file__).resolve().parents[1]
OUTPUT_ROOT = PACKAGE_ROOT / "outputs"
METHOD_ALIASES = {"DM": "DM-ESPCA", "AWGE": "AWGE-ESPCA", "GSN": "GSN-SPCA"}
SUMMARY_METHOD_ALIASES = {"ESPCA": "ESPCA", "DM": "DM-ESPCA", "AWGE": "AWGE-ESPCA", "GSN": "GSN-SPCA", "GSN-SPCA": "GSN-SPCA"}
METHOD_TO_SUMMARY = {"ESPCA": "ESPCA", "DM-ESPCA": "DM", "AWGE-ESPCA": "AWGE", "GSN-SPCA": "GSN"}
CLASS_LABELS = ["PC1", "PC2", "NOISE-LOW", "NOISE-HIGH"]
CLASS_COLORS = {
    "PC1": "#2563eb",
    "PC2": "#059669",
    "NOISE-LOW": "#b7791f",
    "NOISE-HIGH": "#dc2626",
    "BACKGROUND": "#8a928d",
}


def load_config(path: Path) -> dict[str, Any]:
    config = json.loads(path.read_text(encoding="utf-8-sig"))
    variables = {
        "source_root": str(Path(config["source_root"])),
        "package_root": str(PACKAGE_ROOT),
    }
    variables["simulation_root"] = config["simulation_root"].format(**variables)
    variables["ready_dir"] = config["ready_dir"].format(**variables)
    variables["original_pcas_dir"] = config["original_pcas_dir"].format(**variables)
    config["_vars"] = variables
    config.setdefault("parameter_sets_csv", str(PACKAGE_ROOT / "configs" / "parameter_sets.csv"))
    for method_cfg in config["methods"].values():
        for key in ["method_dir", "source_input_dir"]:
            method_cfg[key] = method_cfg[key].format(**variables)
    return config


def select_methods(config: dict[str, Any], methods_arg: str) -> list[str]:
    names = list(config["methods"].keys())
    if methods_arg.lower() == "all":
        return names
    requested = [METHOD_ALIASES.get(item.strip(), item.strip()) for item in methods_arg.split(",") if item.strip()]
    missing = [name for name in requested if name not in config["methods"]]
    if missing:
        raise SystemExit(f"Unknown method(s): {', '.join(missing)}")
    return requested


def read_csv_dicts(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def select_parameter_sets(config: dict[str, Any], sets_arg: str) -> list[dict[str, str]]:
    csv_path = Path(config.get("parameter_sets_csv", PACKAGE_ROOT / "configs" / "parameter_sets.csv"))
    if not csv_path.exists() or sets_arg.lower() in {"", "default", "none"}:
        return []
    rows = read_csv_dicts(csv_path)
    if sets_arg.lower() == "all":
        return rows
    wanted = [item.strip() for item in sets_arg.split(",") if item.strip()]
    by_id = {row["config_id"]: row for row in rows}
    missing = [item for item in wanted if item not in by_id]
    if missing:
        raise SystemExit(f"Unknown parameter_set(s): {', '.join(missing)}")
    return [by_id[item] for item in wanted]


def make_run_dir(run_id: str | None) -> Path:
    if not run_id:
        run_id = "run_" + datetime.now().strftime("%Y%m%d_%H%M%S")
    run_dir = OUTPUT_ROOT / run_id
    if run_dir.exists():
        suffix = datetime.now().strftime("_%f")
        run_dir = OUTPUT_ROOT / f"{run_id}{suffix}"
    (run_dir / "methods").mkdir(parents=True, exist_ok=True)
    (run_dir / "logs").mkdir(parents=True, exist_ok=True)
    (run_dir / "tables").mkdir(parents=True, exist_ok=True)
    return run_dir


def copy_tree_files(src_dir: Path, dest_dir: Path) -> list[str]:
    copied = []
    dest_dir.mkdir(parents=True, exist_ok=True)
    for src in src_dir.iterdir():
        if src.is_file():
            shutil.copy2(src, dest_dir / src.name)
            copied.append(src.name)
    return copied


def copy_method_code(method_cfg: dict[str, Any], dest_dir: Path) -> list[str]:
    return copy_tree_files(Path(method_cfg["method_dir"]), dest_dir)


def patch_gsn_python(method: str, work_dir: Path, python_path: str | None) -> None:
    if method != "GSN-SPCA":
        return
    fun_path = work_dir / "fun_GSN-SPCA.R"
    if not fun_path.exists():
        return
    text = fun_path.read_text(encoding="utf-8", errors="replace")
    if python_path:
        normalized_python = python_path.replace("\\", "/")
        replacement = f'use_python("{normalized_python}", required = TRUE)'
    else:
        replacement = 'use_python(Sys.getenv("RETICULATE_PYTHON"), required = FALSE)'
    lines = []
    replaced = False
    for line in text.splitlines():
        if line.strip().startswith("use_python("):
            lines.append(replacement)
            replaced = True
        else:
            lines.append(line)
    if not replaced:
        lines.insert(1, replacement)
    fun_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def v_matrix_dimensions(path: Path) -> tuple[int | None, int | None]:
    if not path.exists():
        return None, None
    lines = [line for line in path.read_text(encoding="utf-8-sig", errors="replace").splitlines() if line.strip()]
    if not lines:
        return 0, 0
    header = lines[0].split("\t")
    data = lines[1:] if header and header[0] in {"V1", "V2", "x"} else lines
    if not data:
        return 0, len(header)
    return len(data), len(data[0].split("\t"))


def copy_default_inputs(method_cfg: dict[str, Any], dest_dir: Path, prepared_inputs: Path | None, method: str) -> tuple[list[str], list[str]]:
    src_dir = prepared_inputs / method if prepared_inputs else Path(method_cfg["source_input_dir"])
    copied, missing = [], []
    for filename in method_cfg.get("input_files", []):
        src = src_dir / filename
        if src.exists():
            shutil.copy2(src, dest_dir / filename)
            copied.append(filename)
        else:
            missing.append(str(src))
    return copied, missing


def copy_parameter_set_inputs(parameter_set: dict[str, str], method: str, work_dir: Path) -> tuple[list[str], list[str]]:
    input_dir = Path(parameter_set["input_dir"])
    copied, missing = [], []
    mapping = {
        "matrix_x.txt": "gene_new_o.txt",
        "metadata.json": "metadata.json",
        "pc_loadings_generated.csv": "pc_loadings_generated.csv",
    }
    if method == "GSN-SPCA":
        standard_groups = input_dir / "groups_gsn_standard_1based.txt"
        mapping["groups_gsn_standard_1based.txt" if standard_groups.exists() else "groups_gsn_1based.txt"] = "result-1_p2.txt"
    else:
        mapping["edges_1based.txt"] = "result-1_p2.txt"
    for src_name, dest_name in mapping.items():
        src = input_dir / src_name
        if src.exists():
            shutil.copy2(src, work_dir / dest_name)
            copied.append(f"{src_name}->{dest_name}")
        else:
            missing.append(str(src))
    return copied, missing


def summary_lookup(config: dict[str, Any]) -> dict[tuple[str, str], dict[str, str]]:
    path_value = config.get("parameter_set_collect_summary")
    if not path_value:
        return {}
    path = Path(path_value)
    if not path.exists():
        return {}
    out = {}
    for row in read_csv_dicts(path):
        method = SUMMARY_METHOD_ALIASES.get(row.get("method", ""), row.get("method", ""))
        out[(row.get("config_id", ""), method)] = row
    return out


def resolve_attempt_dir(path_text: str) -> Path:
    path = Path(path_text)
    if path.exists():
        return path
    # Fallback: if path_text contains a known marker, look under package root
    text = path_text.replace("/", "\\")
    markers = ["simulation_runs_", "simulation_gsn_baseline_audit_"]
    for marker in markers:
        idx = text.find(marker)
        if idx >= 0:
            candidate = Path(PACKAGE_ROOT).resolve() / text[idx:]
            if candidate.exists():
                return candidate
    return path


def collect_from_attempt(parameter_set: dict[str, str], method: str, summary_row: dict[str, str], dest_dir: Path) -> dict[str, Any]:
    src_dir = resolve_attempt_dir(summary_row.get("attempt_dir", ""))
    copied, missing = [], []
    for filename in ["U_matrix.txt", "D_matrix.txt", "V_matrix.txt", "metadata.json", "attempt_metadata.json", "power.txt"]:
        src = src_dir / filename
        if src.exists():
            shutil.copy2(src, dest_dir / filename)
            copied.append(filename)
        elif filename in {"U_matrix.txt", "D_matrix.txt", "V_matrix.txt"}:
            missing.append(str(src))
    rows, cols = v_matrix_dimensions(dest_dir / "V_matrix.txt")
    return {
        "parameter_set": parameter_set.get("config_id", "default"),
        "family": parameter_set.get("family", ""),
        "method": method,
        "kg": summary_row.get("kg", ""),
        "returncode": 0 if not missing else 1,
        "elapsed_seconds": summary_row.get("elapsed_seconds", 0),
        "stdout": "",
        "stderr": "",
        "work_dir": str(dest_dir),
        "source_dir": str(src_dir),
        "copied_outputs": copied,
        "missing_outputs": missing,
        "V_rows": rows,
        "V_cols": cols,
        "target_hit_count": summary_row.get("target_hit_count", ""),
        "target_hit_rate": summary_row.get("target_hit_rate", ""),
        "noise_selected_count": summary_row.get("noise_selected_count", ""),
        "selected_total": summary_row.get("selected_total", ""),
    }


def collect_default_outputs(method: str, method_cfg: dict[str, Any], dest_dir: Path) -> dict[str, Any]:
    src_dir = Path(method_cfg["source_input_dir"])
    copied, missing = [], []
    for filename in ["U_matrix.txt", "D_matrix.txt", "V_matrix.txt"]:
        src = src_dir / filename
        if src.exists():
            shutil.copy2(src, dest_dir / filename)
            copied.append(filename)
        else:
            missing.append(str(src))
    rows, cols = v_matrix_dimensions(dest_dir / "V_matrix.txt")
    return {
        "parameter_set": "default",
        "family": "default",
        "method": method,
        "returncode": 0 if not missing else 1,
        "elapsed_seconds": 0,
        "stdout": "",
        "stderr": "",
        "work_dir": str(dest_dir),
        "source_dir": str(src_dir),
        "copied_outputs": copied,
        "missing_outputs": missing,
        "V_rows": rows,
        "V_cols": cols,
    }


def numeric_param(value: Any, fallback: float | int) -> Any:
    if value == "auto_from_summary":
        return fallback
    return value


def run_one_method(config: dict[str, Any], method: str, method_cfg: dict[str, Any], work_dir: Path, log_dir: Path, kg_override: int | None) -> dict[str, Any]:
    params = dict(method_cfg["params"])
    if kg_override is not None:
        params["k_group"] = kg_override
    params["k_group"] = numeric_param(params.get("k_group"), method_cfg.get("hardcoded_edge_n") or 1)
    rscript = Path(config["rscript"])
    if not rscript.exists():
        raise RuntimeError(f"Rscript not found: {rscript}")
    runner = PACKAGE_ROOT / "scripts" / "run_method.R"
    cmd = [
        str(rscript), str(runner),
        "--method", method,
        "--work-dir", str(work_dir),
        "--k", str(params["k"]),
        "--k-group", str(params["k_group"]),
        "--we", str(params["we"]),
        "--t", str(params["t"]),
        "--niter", str(params["niter"]),
        "--err", str(params["err"]),
        "--num-init", str(params["num_init"]),
        "--edge-mode", method_cfg.get("edge_mode", "nrow"),
    ]
    if method_cfg.get("hardcoded_edge_n") is not None:
        cmd.extend(["--hardcoded-edge-n", str(method_cfg["hardcoded_edge_n"])])
    safe_name = work_dir.parent.name + "_" + method.replace("/", "_")
    stdout_path = log_dir / f"{safe_name}_stdout.txt"
    stderr_path = log_dir / f"{safe_name}_stderr.txt"
    start = datetime.now()
    env = dict(os.environ)
    if config.get("reticulate_python"):
        env["RETICULATE_PYTHON"] = config["reticulate_python"]
    with stdout_path.open("w", encoding="utf-8", errors="replace") as stdout, stderr_path.open("w", encoding="utf-8", errors="replace") as stderr:
        proc = subprocess.run(cmd, cwd=work_dir, stdout=stdout, stderr=stderr, text=True, env=env)
    elapsed = (datetime.now() - start).total_seconds()
    rows, cols = v_matrix_dimensions(work_dir / "V_matrix.txt")
    return {
        "method": method,
        "kg": str(params["k_group"]),
        "returncode": proc.returncode,
        "elapsed_seconds": elapsed,
        "stdout": str(stdout_path),
        "stderr": str(stderr_path),
        "work_dir": str(work_dir),
        "V_rows": rows,
        "V_cols": cols,
    }


def read_clusters(path: Path) -> list[list[int]]:
    clusters = []
    if not path.exists():
        return clusters
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        values = [int(value) + 1 for value in line.split() if value.strip()]
        if values:
            clusters.append(values)
    return clusters


def read_fixed(path: Path) -> tuple[list[int], list[int]]:
    if not path.exists():
        return [], []
    values = [int(value) + 1 for value in path.read_text(encoding="utf-8", errors="replace").split()]
    return [value for value in values if value <= 25], [value for value in values if value > 25]


def read_v_matrix(path: Path) -> list[dict[str, float]]:
    with path.open("r", newline="", encoding="utf-8-sig", errors="replace") as handle:
        return [{"V1": float(row["V1"]), "V2": float(row["V2"])} for row in csv.DictReader(handle, delimiter="\t")]


def chunk(values: list[int], size: int) -> list[list[int]]:
    return [values[index:index + size] for index in range(0, len(values), size)]


def flatten_ints(values: Any) -> list[int]:
    out: list[int] = []
    if values is None:
        return out
    if isinstance(values, list):
        for item in values:
            out.extend(flatten_ints(item))
    elif isinstance(values, tuple):
        for item in values:
            out.extend(flatten_ints(item))
    elif isinstance(values, (int, float)):
        out.append(int(values))
    elif isinstance(values, str) and values.strip():
        out.append(int(float(values)))
    return out


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8-sig"))


def read_pc_loading_classes(input_dir: Path, n_genes: int) -> tuple[set[int], set[int]]:
    pc_path = input_dir / "pc_loadings_generated.csv"
    pc1: set[int] = set()
    pc2: set[int] = set()
    if not pc_path.exists():
        return pc1, pc2
    with pc_path.open("r", newline="", encoding="utf-8-sig", errors="replace") as handle:
        for index, row in enumerate(csv.DictReader(handle), start=1):
            if n_genes and index > n_genes:
                break
            u1 = float(row.get("u1") or row.get("V1") or 0)
            u2 = float(row.get("u2") or row.get("V2") or 0)
            has_u1 = abs(u1) > 1e-12
            has_u2 = abs(u2) > 1e-12
            if has_u1 and (not has_u2 or abs(u1) >= abs(u2)):
                pc1.add(index)
            elif has_u2:
                pc2.add(index)
    return pc1, pc2


def class_definition_from_parameter_set(parameter_set: dict[str, str]) -> dict[str, Any]:
    input_dir = Path(parameter_set.get("input_dir", ""))
    metadata_value = parameter_set.get("metadata_path", "")
    metadata_path = Path(metadata_value) if metadata_value else input_dir / "metadata.json"
    metadata = load_json(metadata_path)
    config = metadata.get("config", {})
    n_genes = int(float(parameter_set.get("n_genes") or config.get("n_genes") or 0))

    clusters = metadata.get("target_clusters_0based", [])
    midpoint = len(clusters) // 2
    pc1 = {value + 1 for value in flatten_ints(clusters[:midpoint])}
    pc2 = {value + 1 for value in flatten_ints(clusters[midpoint:])}

    noise_low = {
        value + 1
        for value in flatten_ints(metadata.get("fixed_cluster_0based", []))
        if value < 50
    }
    noise_high = {
        value + 1
        for value in flatten_ints(metadata.get("added_noise_clusters_0based", []))
        if value >= 50
    }
    return {
        "n_genes": n_genes,
        "pc1": pc1,
        "pc2": pc2,
        "noise_low": noise_low,
        "noise_high": noise_high,
        "source": str(metadata_path),
    }


def default_class_definition(config: dict[str, Any]) -> dict[str, Any]:
    ready_dir = Path(config["_vars"]["ready_dir"])
    fixed_group1, fixed_group2 = read_fixed(ready_dir / "fixed_cluster_merged.txt")
    return {
        "n_genes": 50,
        "pc1": set(flatten_ints(read_clusters(ready_dir / "clusters_group1_remaining.txt"))),
        "pc2": set(flatten_ints(read_clusters(ready_dir / "clusters_group2_remaining.txt"))),
        "noise_low": set(fixed_group1 + fixed_group2),
        "noise_high": set(),
        "source": str(ready_dir),
    }


def load_class_definitions(config: dict[str, Any]) -> dict[str, dict[str, Any]]:
    definitions = {"default": default_class_definition(config)}
    csv_path = Path(config.get("parameter_sets_csv", PACKAGE_ROOT / "configs" / "parameter_sets.csv"))
    if csv_path.exists():
        for row in read_csv_dicts(csv_path):
            definitions[row["config_id"]] = class_definition_from_parameter_set(row)
    return definitions


def classify_gene(gene_index_1based: int, definition: dict[str, Any]) -> str:
    if gene_index_1based in definition.get("pc1", set()):
        return "PC1"
    if gene_index_1based in definition.get("pc2", set()):
        return "PC2"
    if gene_index_1based in definition.get("noise_low", set()):
        return "NOISE-LOW"
    if gene_index_1based in definition.get("noise_high", set()):
        return "NOISE-HIGH"
    return "BACKGROUND"


def class_counts(definition: dict[str, Any], n_rows: int) -> dict[str, int]:
    counts = {label: 0 for label in CLASS_LABELS}
    counts["BACKGROUND"] = 0
    for gene_index in range(1, n_rows + 1):
        counts[classify_gene(gene_index, definition)] += 1
    return counts


def emit_panel(writer: csv.DictWriter, parameter_set: str, method: str, matrix: list[dict[str, float]], panel: str, rows: list[list[int]], component: str) -> None:
    for panel_row, source_rows in enumerate(rows, start=1):
        for panel_col, source_row in enumerate(source_rows, start=1):
            if source_row - 1 >= len(matrix):
                continue
            writer.writerow({
                "parameter_set": parameter_set,
                "method": method,
                "panel": panel,
                "panel_row": panel_row,
                "panel_col": panel_col,
                "source_row_1based": source_row,
                "source_index_0based": source_row - 1,
                "source_component": component,
                "value": matrix[source_row - 1][component],
            })


def safe_filename(value: str) -> str:
    return "".join(char if char.isalnum() or char in {"-", "_", "."} else "_" for char in value)


def opacity_for_value(value: float, max_abs: float) -> float:
    if max_abs <= 0:
        return 0.18
    return 0.18 + 0.78 * min(abs(value) / max_abs, 1.0)


def make_bins(values: list[int], max_bins: int) -> list[list[int]]:
    if len(values) <= max_bins:
        return [[value] for value in values]
    bins: list[list[int]] = []
    total = len(values)
    for index in range(max_bins):
        start = index * total // max_bins
        end = (index + 1) * total // max_bins
        part = values[start:end]
        if part:
            bins.append(part)
    return bins


def write_loading_svg(run_dir: Path, result: dict[str, Any], matrix: list[dict[str, float]], definition: dict[str, Any]) -> dict[str, Any]:
    figures_dir = run_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    parameter_set = str(result.get("parameter_set", "default"))
    method = str(result.get("method", "method"))
    out_path = figures_dir / f"{safe_filename(parameter_set)}__{safe_filename(method)}__classes.svg"
    n_rows = len(matrix)
    max_abs = max([abs(row["V1"]) for row in matrix] + [abs(row["V2"]) for row in matrix] + [0.0])
    counts = class_counts(definition, n_rows)
    outside_count = max(n_rows - 50, 0)
    height = 500 if outside_count else 390
    width = 1280
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect x="0" y="0" width="1280" height="100%" fill="#fbfbf8"/>',
        '<style>text{font-family:Segoe UI,Arial,sans-serif;fill:#202926}.title{font-size:22px;font-weight:700}.sub{font-size:13px;fill:#4b5852}.axis{font-size:11px;fill:#52615a}.label{font-size:13px;font-weight:600}.small{font-size:11px;fill:#425049}</style>',
        f'<text x="28" y="34" class="title">{html.escape(parameter_set)} / {html.escape(method)}</text>',
        f'<text x="28" y="58" class="sub">Class labels: PC1 and PC2 are metadata target clusters; NOISE-LOW is the first-50 fixed-noise set; NOISE-HIGH is outside the first 50 genes.</text>',
    ]
    legend_x = 28
    for label in CLASS_LABELS:
        color = CLASS_COLORS[label]
        parts.append(f'<rect x="{legend_x}" y="78" width="14" height="14" rx="2" fill="{color}"/>')
        parts.append(f'<text x="{legend_x + 20}" y="90" class="axis">{label} ({counts.get(label, 0)})</text>')
        legend_x += 150
    parts.append(f'<text x="{legend_x}" y="90" class="axis">BACKGROUND ({counts.get("BACKGROUND", 0)})</text>')

    left = 96
    top = 132
    cell_w = 17
    cell_h = 24
    gap = 2
    inside = list(range(1, min(50, n_rows) + 1))
    parts.append(f'<text x="28" y="{top - 18}" class="label">Inside first 50 genes</text>')
    for row_index, component in enumerate(["V1", "V2"]):
        y = top + row_index * (cell_h + 8)
        parts.append(f'<text x="58" y="{y + 16}" class="axis">{component}</text>')
        for offset, gene_index in enumerate(inside):
            x = left + offset * (cell_w + gap)
            label = classify_gene(gene_index, definition)
            value = matrix[gene_index - 1][component]
            opacity = opacity_for_value(value, max_abs)
            color = CLASS_COLORS.get(label, CLASS_COLORS["BACKGROUND"])
            parts.append(f'<rect x="{x}" y="{y}" width="{cell_w}" height="{cell_h}" rx="2" fill="{color}" fill-opacity="{opacity:.3f}" stroke="#ffffff" stroke-width="0.6"><title>Gene {gene_index}; {label}; {component}={value:.6g}</title></rect>')
    for gene_index in [1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50]:
        if gene_index <= len(inside):
            x = left + (gene_index - 1) * (cell_w + gap)
            parts.append(f'<text x="{x}" y="{top + 2 * (cell_h + 8) + 12}" class="axis">{gene_index}</text>')

    outside = list(range(51, n_rows + 1))
    if outside:
        outside_top = top + 128
        bins = make_bins(outside, 120)
        bin_w = 8
        parts.append(f'<text x="28" y="{outside_top - 20}" class="label">Outside first 50 genes</text>')
        parts.append(f'<text x="220" y="{outside_top - 20}" class="sub">{outside_count} genes displayed as {len(bins)} ordered bins; color is NOISE-HIGH when defined by metadata.</text>')
        for row_index, component in enumerate(["V1", "V2"]):
            y = outside_top + row_index * (cell_h + 8)
            parts.append(f'<text x="58" y="{y + 16}" class="axis">{component}</text>')
            for bin_index, genes in enumerate(bins):
                values = [matrix[gene - 1][component] for gene in genes if gene - 1 < len(matrix)]
                value = max(values, key=lambda item: abs(item)) if values else 0.0
                label = "NOISE-HIGH" if any(classify_gene(gene, definition) == "NOISE-HIGH" for gene in genes) else "BACKGROUND"
                color = CLASS_COLORS.get(label, CLASS_COLORS["BACKGROUND"])
                opacity = opacity_for_value(value, max_abs)
                x = left + bin_index * (bin_w + 1)
                if len(genes) == 1:
                    title = f"Gene {genes[0]}; {label}; {component}={value:.6g}"
                else:
                    title = f"Genes {genes[0]}-{genes[-1]}; {label}; max abs {component}={value:.6g}"
                parts.append(f'<rect x="{x}" y="{y}" width="{bin_w}" height="{cell_h}" rx="1" fill="{color}" fill-opacity="{opacity:.3f}" stroke="#ffffff" stroke-width="0.35"><title>{html.escape(title)}</title></rect>')
        parts.append(f'<text x="{left}" y="{outside_top + 2 * (cell_h + 8) + 12}" class="axis">51</text>')
        parts.append(f'<text x="{left + (len(bins) - 1) * (bin_w + 1) - 10}" y="{outside_top + 2 * (cell_h + 8) + 12}" class="axis">{n_rows}</text>')
    else:
        parts.append(f'<text x="28" y="{top + 126}" class="sub">No genes outside the first 50 are present in this V matrix.</text>')

    parts.append(f'<text x="28" y="{height - 24}" class="small">Opacity is proportional to absolute loading magnitude. Full numeric values are available in tables/classified_loadings.csv.</text>')
    parts.append("</svg>")
    out_path.write_text("\n".join(parts), encoding="utf-8")
    return {
        "parameter_set": parameter_set,
        "method": method,
        "file": str(out_path),
        "n_genes": n_rows,
        "inside_50": min(50, n_rows),
        "outside_50": outside_count,
        "counts": counts,
    }


def selected_count_for_class(matrix: list[dict[str, float]], definition: dict[str, Any], label: str) -> int:
    count = 0
    for gene_index, values in enumerate(matrix, start=1):
        if classify_gene(gene_index, definition) != label:
            continue
        if label == "PC1":
            selected = abs(values["V1"]) > 1e-12
        elif label == "PC2":
            selected = abs(values["V2"]) > 1e-12
        else:
            selected = abs(values["V1"]) > 1e-12 or abs(values["V2"]) > 1e-12
        if selected:
            count += 1
    return count


def stacked_bar_rows_for_result(result: dict[str, Any], matrix: list[dict[str, float]], definition: dict[str, Any]) -> list[dict[str, Any]]:
    rows = []
    for label in CLASS_LABELS:
        rows.append({
            "parameter_set": result.get("parameter_set", "default"),
            "method": result.get("method", ""),
            "class_label": label,
            "selected_count": selected_count_for_class(matrix, definition, label),
        })
    return rows


def write_stacked_bar_svg(run_dir: Path, parameter_set: str, rows: list[dict[str, Any]]) -> dict[str, Any]:
    figures_dir = run_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    out_path = figures_dir / f"{safe_filename(parameter_set)}__method_class_counts.svg"
    methods = []
    for row in rows:
        method = str(row["method"])
        if method not in methods:
            methods.append(method)
    totals = {
        method: sum(int(row["selected_count"]) for row in rows if row["method"] == method)
        for method in methods
    }
    max_total = max(totals.values(), default=1)
    width = max(780, 190 + 135 * max(len(methods), 1))
    height = 430
    plot_left = 90
    plot_top = 82
    plot_height = 250
    bar_w = 56
    gap = 76
    scale = plot_height / max(max_total, 1)
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect x="0" y="0" width="100%" height="100%" fill="#fbfbf8"/>',
        '<style>text{font-family:Segoe UI,Arial,sans-serif;fill:#202926}.title{font-size:22px;font-weight:700}.sub{font-size:13px;fill:#4b5852}.axis{font-size:11px;fill:#52615a}.tick{stroke:#d7ddd8;stroke-width:1}.label{font-size:12px;font-weight:600}</style>',
        f'<text x="28" y="34" class="title">{html.escape(parameter_set)}: selected genes by method</text>',
        '<text x="28" y="58" class="sub">Stacked counts: PC1 uses nonzero V1; PC2 uses nonzero V2; NOISE-LOW and NOISE-HIGH use nonzero V1 or V2.</text>',
    ]
    legend_x = 28
    for label in CLASS_LABELS:
        parts.append(f'<rect x="{legend_x}" y="374" width="14" height="14" rx="2" fill="{CLASS_COLORS[label]}"/>')
        parts.append(f'<text x="{legend_x + 20}" y="386" class="axis">{label}</text>')
        legend_x += 150
    tick_step = 1 if max_total <= 8 else max(1, round(max_total / 5))
    tick = 0
    while tick <= max_total:
        y = plot_top + plot_height - tick * scale
        parts.append(f'<line x1="{plot_left - 8}" y1="{y:.1f}" x2="{width - 32}" y2="{y:.1f}" class="tick"/>')
        parts.append(f'<text x="{plot_left - 36}" y="{y + 4:.1f}" class="axis">{tick}</text>')
        tick += tick_step
    if (tick - tick_step) != max_total:
        y = plot_top
        parts.append(f'<line x1="{plot_left - 8}" y1="{y:.1f}" x2="{width - 32}" y2="{y:.1f}" class="tick"/>')
        parts.append(f'<text x="{plot_left - 36}" y="{y + 4:.1f}" class="axis">{max_total}</text>')
    parts.append(f'<line x1="{plot_left - 8}" y1="{plot_top + plot_height}" x2="{width - 32}" y2="{plot_top + plot_height}" stroke="#8a928d" stroke-width="1"/>')
    for method_index, method in enumerate(methods):
        x = plot_left + method_index * (bar_w + gap)
        y_cursor = plot_top + plot_height
        for label in CLASS_LABELS:
            value = sum(int(row["selected_count"]) for row in rows if row["method"] == method and row["class_label"] == label)
            segment_h = value * scale
            y_cursor -= segment_h
            if value > 0:
                parts.append(f'<rect x="{x}" y="{y_cursor:.1f}" width="{bar_w}" height="{segment_h:.1f}" fill="{CLASS_COLORS[label]}" stroke="#ffffff" stroke-width="0.8"><title>{html.escape(method)}; {label}; selected_count={value}</title></rect>')
                if segment_h >= 16:
                    parts.append(f'<text x="{x + bar_w / 2:.1f}" y="{y_cursor + segment_h / 2 + 4:.1f}" text-anchor="middle" class="axis" fill="#ffffff">{value}</text>')
        total = totals.get(method, 0)
        parts.append(f'<text x="{x + bar_w / 2:.1f}" y="{plot_top + plot_height + 20}" text-anchor="middle" class="label">{html.escape(method)}</text>')
        parts.append(f'<text x="{x + bar_w / 2:.1f}" y="{max(y_cursor - 8, 72):.1f}" text-anchor="middle" class="axis">total {total}</text>')
    parts.append("</svg>")
    out_path.write_text("\n".join(parts), encoding="utf-8")
    return {"parameter_set": parameter_set, "file": str(out_path), "methods": methods, "max_total": max_total}


def class_value(rows: list[dict[str, Any]], parameter_set: str, method: str, label: str) -> int:
    for row in rows:
        if row["parameter_set"] == parameter_set and row["method"] == method and row["class_label"] == label:
            return int(row["selected_count"])
    return 0


def ordered_unique(values: list[str]) -> list[str]:
    out = []
    for value in values:
        if value not in out:
            out.append(value)
    return out


def write_progress_table_svg(run_dir: Path, stacked_rows: list[dict[str, Any]]) -> dict[str, Any]:
    figures_dir = run_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    out_path = figures_dir / "overview_10x4_progress_table.svg"
    parameter_sets = ordered_unique([str(row["parameter_set"]) for row in stacked_rows])
    methods = ordered_unique([str(row["method"]) for row in stacked_rows])
    totals = {
        (ps, method): sum(class_value(stacked_rows, ps, method, label) for label in CLASS_LABELS)
        for ps in parameter_sets
        for method in methods
    }
    max_total = max(totals.values(), default=1)
    left = 260
    top = 138
    row_h = 34
    col_w = 245
    bar_w = 185
    bar_h = 16
    width = left + col_w * max(len(methods), 1) + 30
    height = top + row_h * max(len(parameter_sets), 1) + 76
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect x="0" y="0" width="100%" height="100%" fill="#fbfbf8"/>',
        '<style>text{font-family:Segoe UI,Arial,sans-serif;fill:#202926}.title{font-size:22px;font-weight:700}.sub{font-size:13px;fill:#4b5852}.head{font-size:12px;font-weight:700}.cell{font-size:11px;fill:#52615a}.rowline{stroke:#e2e6e2;stroke-width:1}</style>',
        '<text x="28" y="34" class="title">10 x 4 Method Progress Table</text>',
        '<text x="28" y="58" class="sub">Each cell is a horizontal stacked bar using selected-gene counts from stacked_bar_counts.csv.</text>',
    ]
    legend_x = 28
    for label in CLASS_LABELS:
        parts.append(f'<rect x="{legend_x}" y="78" width="13" height="13" rx="2" fill="{CLASS_COLORS[label]}"/>')
        parts.append(f'<text x="{legend_x + 18}" y="89" class="cell">{label}</text>')
        legend_x += 145
    for method_index, method in enumerate(methods):
        x = left + method_index * col_w
        parts.append(f'<text x="{x}" y="{top - 18}" class="head">{html.escape(method)}</text>')
    for row_index, ps in enumerate(parameter_sets):
        y = top + row_index * row_h
        parts.append(f'<line x1="24" y1="{y + 24}" x2="{width - 24}" y2="{y + 24}" class="rowline"/>')
        parts.append(f'<text x="28" y="{y + 14}" class="cell">{html.escape(ps)}</text>')
        for method_index, method in enumerate(methods):
            x = left + method_index * col_w
            total = totals.get((ps, method), 0)
            denom = max(max_total, 1)
            cursor = x
            for label in CLASS_LABELS:
                value = class_value(stacked_rows, ps, method, label)
                seg_w = bar_w * value / denom
                if seg_w > 0:
                    parts.append(f'<rect x="{cursor:.1f}" y="{y}" width="{seg_w:.1f}" height="{bar_h}" fill="{CLASS_COLORS[label]}"><title>{html.escape(ps)}; {html.escape(method)}; {label}; selected_count={value}; total={total}</title></rect>')
                cursor += seg_w
            parts.append(f'<rect x="{x}" y="{y}" width="{bar_w}" height="{bar_h}" fill="none" stroke="#76827b" stroke-width="0.6"/>')
            parts.append(f'<text x="{x + bar_w + 7}" y="{y + 12}" class="cell">{total}</text>')
    parts.append(f'<text x="28" y="{height - 24}" class="sub">Bar length is scaled to the largest method total in this run: {max_total}.</text>')
    parts.append("</svg>")
    out_path.write_text("\n".join(parts), encoding="utf-8")
    return {"file": str(out_path), "parameter_sets": len(parameter_sets), "methods": len(methods), "max_total": max_total}


def write_method_parameter_svg(run_dir: Path, method: str, stacked_rows: list[dict[str, Any]]) -> dict[str, Any]:
    figures_dir = run_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    out_path = figures_dir / f"method_{safe_filename(method)}_by_parameter_set.svg"
    parameter_sets = ordered_unique([str(row["parameter_set"]) for row in stacked_rows])
    totals = {
        ps: sum(class_value(stacked_rows, ps, method, label) for label in CLASS_LABELS)
        for ps in parameter_sets
    }
    max_total = max(totals.values(), default=1)
    left = 260
    top = 82
    row_h = 38
    bar_w = 620
    bar_h = 20
    width = 970
    height = top + row_h * max(len(parameter_sets), 1) + 72
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect x="0" y="0" width="100%" height="100%" fill="#fbfbf8"/>',
        '<style>text{font-family:Segoe UI,Arial,sans-serif;fill:#202926}.title{font-size:22px;font-weight:700}.sub{font-size:13px;fill:#4b5852}.cell{font-size:11px;fill:#52615a}.rowline{stroke:#e2e6e2;stroke-width:1}</style>',
        f'<text x="28" y="34" class="title">{html.escape(method)} by Parameter Set</text>',
        '<text x="28" y="58" class="sub">Horizontal stacked bars compare selected PC/noise counts across the 10 parameter sets.</text>',
    ]
    legend_x = 28
    for label in CLASS_LABELS:
        parts.append(f'<rect x="{legend_x}" y="{height - 34}" width="13" height="13" rx="2" fill="{CLASS_COLORS[label]}"/>')
        parts.append(f'<text x="{legend_x + 18}" y="{height - 23}" class="cell">{label}</text>')
        legend_x += 145
    for row_index, ps in enumerate(parameter_sets):
        y = top + row_index * row_h
        total = totals.get(ps, 0)
        cursor = left
        parts.append(f'<line x1="24" y1="{y + 28}" x2="{width - 24}" y2="{y + 28}" class="rowline"/>')
        parts.append(f'<text x="28" y="{y + 15}" class="cell">{html.escape(ps)}</text>')
        for label in CLASS_LABELS:
            value = class_value(stacked_rows, ps, method, label)
            seg_w = bar_w * value / max(max_total, 1)
            if seg_w > 0:
                parts.append(f'<rect x="{cursor:.1f}" y="{y}" width="{seg_w:.1f}" height="{bar_h}" fill="{CLASS_COLORS[label]}"><title>{html.escape(method)}; {html.escape(ps)}; {label}; selected_count={value}; total={total}</title></rect>')
            cursor += seg_w
        parts.append(f'<rect x="{left}" y="{y}" width="{bar_w}" height="{bar_h}" fill="none" stroke="#76827b" stroke-width="0.6"/>')
        parts.append(f'<text x="{left + bar_w + 12}" y="{y + 15}" class="cell">total {total}</text>')
    parts.append("</svg>")
    out_path.write_text("\n".join(parts), encoding="utf-8")
    return {"method": method, "file": str(out_path), "parameter_sets": len(parameter_sets), "max_total": max_total}


def write_method_parameter_figures(run_dir: Path, stacked_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    methods = ordered_unique([str(row["method"]) for row in stacked_rows])
    return [write_method_parameter_svg(run_dir, method, stacked_rows) for method in methods]


def write_stacked_bar_figures(run_dir: Path, stacked_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    figures = []
    parameter_sets = []
    for row in stacked_rows:
        ps = str(row["parameter_set"])
        if ps not in parameter_sets:
            parameter_sets.append(ps)
    for ps in parameter_sets:
        rows = [row for row in stacked_rows if row["parameter_set"] == ps]
        figures.append(write_stacked_bar_svg(run_dir, ps, rows))
    return figures


def build_heatmap_values(config: dict[str, Any], run_dir: Path, results: list[dict[str, Any]]) -> dict[str, Any]:
    definitions = load_class_definitions(config)
    out_path = run_dir / "tables" / "default_heatmap_values.csv"
    classified_path = run_dir / "tables" / "classified_loadings.csv"
    stacked_path = run_dir / "tables" / "stacked_bar_counts.csv"
    heatmap_fields = ["parameter_set", "method", "panel", "panel_row", "panel_col", "source_row_1based", "source_index_0based", "source_component", "value"]
    classified_fields = ["parameter_set", "method", "gene_index_1based", "gene_index_0based", "region", "class_label", "V1", "V2", "abs_max"]
    stacked_fields = ["parameter_set", "method", "class_label", "selected_count"]
    written, missing, class_rows = 0, [], 0
    figures = []
    stacked_rows: list[dict[str, Any]] = []
    with out_path.open("w", newline="", encoding="utf-8-sig") as heatmap_handle, classified_path.open("w", newline="", encoding="utf-8-sig") as class_handle, stacked_path.open("w", newline="", encoding="utf-8-sig") as stacked_handle:
        heatmap_writer = csv.DictWriter(heatmap_handle, fieldnames=heatmap_fields)
        class_writer = csv.DictWriter(class_handle, fieldnames=classified_fields)
        stacked_writer = csv.DictWriter(stacked_handle, fieldnames=stacked_fields)
        heatmap_writer.writeheader()
        class_writer.writeheader()
        stacked_writer.writeheader()
        for result in results:
            v_path = Path(result["work_dir"]) / "V_matrix.txt"
            if not v_path.exists():
                missing.append(str(v_path))
                continue
            matrix = read_v_matrix(v_path)
            ps = result.get("parameter_set", "default")
            method = result["method"]
            definition = definitions.get(ps, definitions["default"])
            for gene_index, values in enumerate(matrix, start=1):
                label = classify_gene(gene_index, definition)
                class_writer.writerow({
                    "parameter_set": ps,
                    "method": method,
                    "gene_index_1based": gene_index,
                    "gene_index_0based": gene_index - 1,
                    "region": "INSIDE-50" if gene_index <= 50 else "OUTSIDE-50",
                    "class_label": label,
                    "V1": values["V1"],
                    "V2": values["V2"],
                    "abs_max": max(abs(values["V1"]), abs(values["V2"])),
                })
                class_rows += 1
            method_stacked_rows = stacked_bar_rows_for_result(result, matrix, definition)
            stacked_rows.extend(method_stacked_rows)
            stacked_writer.writerows(method_stacked_rows)
            emit_panel(heatmap_writer, ps, method, matrix, "PC1", chunk(sorted(definition.get("pc1", [])), 5), "V1")
            emit_panel(heatmap_writer, ps, method, matrix, "PC2", chunk(sorted(definition.get("pc2", [])), 5), "V2")
            emit_panel(heatmap_writer, ps, method, matrix, "NOISE-LOW-V1", chunk(sorted(definition.get("noise_low", [])), 5), "V1")
            emit_panel(heatmap_writer, ps, method, matrix, "NOISE-LOW-V2", chunk(sorted(definition.get("noise_low", [])), 5), "V2")
            figures.append(write_loading_svg(run_dir, result, matrix, definition))
            written += 1
    stacked_figures = write_stacked_bar_figures(run_dir, stacked_rows)
    progress_table_figure = write_progress_table_svg(run_dir, stacked_rows)
    method_parameter_figures = write_method_parameter_figures(run_dir, stacked_rows)
    return {
        "path": str(out_path),
        "classified_path": str(classified_path),
        "stacked_path": str(stacked_path),
        "figures_dir": str(run_dir / "figures"),
        "figures": figures,
        "stacked_figures": stacked_figures,
        "progress_table_figure": progress_table_figure,
        "method_parameter_figures": method_parameter_figures,
        "methods_written": written,
        "classified_rows": class_rows,
        "missing": missing,
        "panel_layout": "PC1 and PC2 are assigned from metadata target clusters. NOISE-LOW is the first-50 fixed-noise set. NOISE-HIGH is the metadata-defined added-noise set outside the first 50 genes.",
    }


def write_method_summary_csv(run_dir: Path, results: list[dict[str, Any]]) -> Path:
    out = run_dir / "tables" / "method_run_summary.csv"
    fieldnames = ["parameter_set", "family", "method", "kg", "returncode", "elapsed_seconds", "V_rows", "V_cols", "target_hit_count", "target_hit_rate", "noise_selected_count", "selected_total", "work_dir", "source_dir", "stdout", "stderr"]
    with out.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(results)
    return out


def write_html(run_dir: Path, mode: str, config: dict[str, Any], parameter_sets: list[dict[str, str]], results: list[dict[str, Any]], heatmap: dict[str, Any]) -> Path:
    rows = []
    for result in results:
        status = "success" if result.get("returncode") == 0 else "failed"
        dims = "missing" if result.get("V_rows") is None else f"{result.get('V_rows')} x {result.get('V_cols')}"
        rows.append("<tr>" +
            f"<td>{html.escape(str(result.get('parameter_set', 'default')))}</td>" +
            f"<td>{html.escape(str(result['method']))}</td>" +
            f"<td>{html.escape(str(result.get('kg', '')))}</td>" +
            f"<td>{html.escape(status)}</td>" +
            f"<td>{html.escape(dims)}</td>" +
            f"<td>{html.escape(str(result.get('target_hit_count', '')))}</td>" +
            f"<td>{html.escape(str(result.get('noise_selected_count', '')))}</td>" +
            f"<td><code>{html.escape(str(result.get('work_dir', '')))}</code></td>" +
            "</tr>")
    ps_rows = []
    for ps in parameter_sets:
        ps_rows.append("<tr>" + "".join(f"<td>{html.escape(str(ps.get(col, '')))}</td>" for col in ["config_id", "family", "n_genes", "d1", "d2", "c", "p_target", "p_noise", "total_edges", "gsn_total_groups"]) + "</tr>")
    seed_rows = "".join(f"<tr><td>{html.escape(str(k))}</td><td><code>{html.escape(str(v))}</code></td></tr>" for k, v in config.get("seeds", {}).items())
    progress_card = ""
    progress_figure = heatmap.get("progress_table_figure", {})
    if progress_figure.get("file"):
        rel = os.path.relpath(progress_figure["file"], run_dir).replace("\\", "/")
        progress_card = (
            '<figure class="fig-card wide-card">' +
            '<figcaption><strong>10 x 4 progress table</strong><br><span>Rows are parameter sets; columns are methods; cells are horizontal stacked bars.</span></figcaption>' +
            f'<a href="{html.escape(rel)}"><img src="{html.escape(rel)}" alt="10 by 4 method progress table"></a>' +
            '</figure>'
        )
    method_cards = []
    for figure in heatmap.get("method_parameter_figures", []):
        rel = os.path.relpath(figure["file"], run_dir).replace("\\", "/")
        method_cards.append(
            '<figure class="fig-card">' +
            f'<figcaption><strong>{html.escape(figure["method"])}</strong><br><span>One method compared across the 10 parameter sets.</span></figcaption>' +
            f'<a href="{html.escape(rel)}"><img src="{html.escape(rel)}" alt="{html.escape(figure["method"])} by parameter set"></a>' +
            '</figure>'
        )
    stacked_cards = []
    for figure in heatmap.get("stacked_figures", []):
        rel = os.path.relpath(figure["file"], run_dir).replace("\\", "/")
        stacked_cards.append(
            '<figure class="fig-card">' +
            f'<figcaption><strong>{html.escape(figure["parameter_set"])}</strong><br><span>Stacked selected-gene counts by method and class.</span></figcaption>' +
            f'<a href="{html.escape(rel)}"><img src="{html.escape(rel)}" alt="{html.escape(figure["parameter_set"])} stacked selected-gene counts"></a>' +
            '</figure>'
        )
    figure_cards = []
    for figure in heatmap.get("figures", []):
        rel = os.path.relpath(figure["file"], run_dir).replace("\\", "/")
        counts = figure.get("counts", {})
        count_text = " | ".join(f"{label}: {counts.get(label, 0)}" for label in CLASS_LABELS)
        figure_cards.append(
            '<figure class="fig-card">' +
            f'<figcaption><strong>{html.escape(figure["parameter_set"])}</strong> / {html.escape(figure["method"])}<br>' +
            f'<span>{html.escape(str(figure["n_genes"]))} genes; inside first 50: {html.escape(str(figure["inside_50"]))}; outside first 50: {html.escape(str(figure["outside_50"]))}</span><br>' +
            f'<span>{html.escape(count_text)}</span></figcaption>' +
            f'<a href="{html.escape(rel)}"><img src="{html.escape(rel)}" alt="{html.escape(figure["parameter_set"])} {html.escape(figure["method"])} class visualization"></a>' +
            '</figure>'
        )
    report = run_dir / "run_result.html"
    generated = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    doc = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>Simulation Engine Run</title>
<style>body{{margin:0;background:#f7f7f4;color:#1f2926;font:14px/1.55 "Segoe UI",Arial,sans-serif}}header{{background:#fff;border-bottom:1px solid #d9ded9;padding:24px 30px}}main{{max-width:1320px;margin:0 auto;padding:22px}}section{{background:#fff;border:1px solid #d9ded9;border-radius:8px;padding:18px;margin:16px 0}}h1{{margin:0 0 8px;font-size:25px}}h2{{margin:0 0 12px;font-size:19px}}table{{width:100%;border-collapse:collapse}}th,td{{border-bottom:1px solid #d9ded9;padding:8px 9px;vertical-align:top;text-align:left}}th{{background:#eef2ef}}code{{font-family:Consolas,"Courier New",monospace;font-size:12px}}a{{color:#177260;text-decoration:none}}.note{{background:#eef5f2;border-left:4px solid #177260;padding:9px 12px}}.fig-grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(420px,1fr));gap:16px}}.fig-card{{margin:0;border:1px solid #d9ded9;border-radius:8px;padding:10px;background:#fff}}.fig-card img{{width:100%;height:auto;border:1px solid #e3e7e3;background:#fbfbf8}}figcaption{{font-size:13px;margin-bottom:8px}}figcaption span{{color:#4b5852}}.wide-card{{grid-column:1/-1}}</style></head>
<body><header><h1>Simulation Engine Run</h1><p>Generated at: {html.escape(generated)}</p><p>Mode: <code>{html.escape(mode)}</code></p><p>Run directory: <code>{html.escape(str(run_dir))}</code></p><div class="note">ESPCA-LEE is treated as the GSN-SPCA experiment line in this workspace. The report labels four gene classes: PC1, PC2, NOISE-LOW, and NOISE-HIGH.</div></header><main>
<section><h2>Parameter Sets</h2><table><thead><tr><th>config_id</th><th>family</th><th>n_genes</th><th>d1</th><th>d2</th><th>c</th><th>p_target</th><th>p_noise</th><th>total_edges</th><th>gsn_total_groups</th></tr></thead><tbody>{''.join(ps_rows) if ps_rows else '<tr><td colspan="10">default</td></tr>'}</tbody></table></section>
<section><h2>Method Results</h2><table><thead><tr><th>Parameter set</th><th>Method</th><th>kg/k.group</th><th>Status</th><th>V_matrix</th><th>Target hit</th><th>Noise selected</th><th>Directory</th></tr></thead><tbody>{''.join(rows)}</tbody></table></section>
<section><h2>Reproducibility Seeds</h2><table><thead><tr><th>Location</th><th>Value</th></tr></thead><tbody>{seed_rows}</tbody></table></section>
<section><h2>10 x 4 Progress Table</h2><p>Rows are the 10 parameter sets, columns are the 4 methods, and every cell is a horizontal stacked progress bar.</p><p>Stacked-bar table: <a href="tables/stacked_bar_counts.csv">tables/stacked_bar_counts.csv</a></p><div class="fig-grid">{progress_card if progress_card else '<p>No progress table was generated.</p>'}</div></section>\n<section><h2>Method-Specific Bar Charts</h2><p>These four charts group results by method, with one horizontal stacked bar for each parameter set.</p><div class="fig-grid">{''.join(method_cards) if method_cards else '<p>No method-specific charts were generated.</p>'}</div></section>\n<section><h2>Stacked Bar Charts</h2><p>Each parameter set has one grouped chart. Each method is one stacked bar: PC1 uses nonzero V1 counts, PC2 uses nonzero V2 counts, and both noise classes use nonzero V1 or V2 counts.</p><div class="fig-grid">{''.join(stacked_cards) if stacked_cards else '<p>No stacked bar charts were generated.</p>'}</div></section>
<section><h2>Visualization Outputs</h2><p>{html.escape(str(heatmap.get('panel_layout', '')))}</p><p>Classified loading table: <a href="tables/classified_loadings.csv">tables/classified_loadings.csv</a></p><p>Panel-value table: <a href="tables/default_heatmap_values.csv">tables/default_heatmap_values.csv</a></p><p>Class figures generated: {html.escape(str(len(heatmap.get('figures', []))))}; stacked charts generated: {html.escape(str(len(heatmap.get('stacked_figures', []))))}; method charts generated: {html.escape(str(len(heatmap.get('method_parameter_figures', []))))}; classified rows: {html.escape(str(heatmap.get('classified_rows', 0)))}</p><div class="fig-grid">{''.join(figure_cards) if figure_cards else '<p>No class figures were generated.</p>'}</div></section>
</main></body></html>"""
    report.write_text(doc, encoding="utf-8")
    return report


def write_manifest(run_dir: Path, mode: str, config_path: Path, methods: list[str], parameter_sets: list[dict[str, str]], results: list[dict[str, Any]], heatmap: dict[str, Any]) -> Path:
    manifest = {"generated_at": datetime.now().isoformat(timespec="seconds"), "package_root": str(PACKAGE_ROOT), "mode": mode, "config": str(config_path.resolve()), "methods": methods, "parameter_sets": parameter_sets, "results": results, "heatmap": heatmap, "terminology": "ESPCA-LEE is treated as the GSN-SPCA experiment line."}
    out = run_dir / "manifest.json"
    out.write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    return out


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default=str(PACKAGE_ROOT / "configs" / "default_config.json"))
    parser.add_argument("--mode", choices=["collect", "run"], default="collect")
    parser.add_argument("--methods", default="all")
    parser.add_argument("--parameter-sets", default="default")
    parser.add_argument("--run-id", default="")
    parser.add_argument("--prepared-inputs", default="")
    args = parser.parse_args()

    config_path = Path(args.config)
    config = load_config(config_path)
    methods = select_methods(config, args.methods)
    parameter_sets = select_parameter_sets(config, args.parameter_sets)
    run_dir = make_run_dir(args.run_id or None)
    prepared_inputs = Path(args.prepared_inputs).resolve() if args.prepared_inputs else None
    summary = summary_lookup(config)

    results: list[dict[str, Any]] = []
    work_items = parameter_sets if parameter_sets else [{"config_id": "default", "family": "default"}]
    for ps in work_items:
        ps_id = ps["config_id"]
        for method in methods:
            method_cfg = config["methods"][method]
            if parameter_sets:
                work_dir = run_dir / ps_id / "methods" / method
            else:
                work_dir = run_dir / "methods" / method
            work_dir.mkdir(parents=True, exist_ok=True)
            copied_code = copy_method_code(method_cfg, work_dir)
            patch_gsn_python(method, work_dir, config.get("reticulate_python"))
            kg_override = None
            if parameter_sets:
                row = summary.get((ps_id, method), {})
                if row.get("kg"):
                    kg_override = int(float(row["kg"]))
            if args.mode == "collect":
                if parameter_sets:
                    row = summary.get((ps_id, method))
                    if row:
                        result = collect_from_attempt(ps, method, row, work_dir)
                    else:
                        result = {"parameter_set": ps_id, "family": ps.get("family", ""), "method": method, "returncode": 3, "work_dir": str(work_dir), "V_rows": None, "V_cols": None, "missing_outputs": ["summary row not found"]}
                else:
                    result = collect_default_outputs(method, method_cfg, work_dir)
                    copied_inputs = []
            else:
                if parameter_sets:
                    copied_inputs, missing_inputs = copy_parameter_set_inputs(ps, method, work_dir)
                else:
                    copied_inputs, missing_inputs = copy_default_inputs(method_cfg, work_dir, prepared_inputs, method)
                if missing_inputs:
                    result = {"parameter_set": ps_id, "family": ps.get("family", ""), "method": method, "returncode": 2, "elapsed_seconds": 0, "stdout": "", "stderr": "", "work_dir": str(work_dir), "missing_inputs": missing_inputs, "V_rows": None, "V_cols": None}
                else:
                    result = run_one_method(config, method, method_cfg, work_dir, run_dir / "logs", kg_override)
                    result["parameter_set"] = ps_id
                    result["family"] = ps.get("family", "")
                    result["copied_inputs"] = copied_inputs
            result["copied_code"] = copied_code
            results.append(result)

    summary_csv = write_method_summary_csv(run_dir, results)
    heatmap = build_heatmap_values(config, run_dir, results)
    report = write_html(run_dir, args.mode, config, parameter_sets, results, heatmap)
    manifest = write_manifest(run_dir, args.mode, config_path, methods, parameter_sets, results, heatmap)

    print(f"Run directory: {run_dir}")
    print(f"Summary CSV: {summary_csv}")
    print(f"Panel values CSV: {heatmap['path']}")
    print(f"Classified loadings CSV: {heatmap.get('classified_path', '')}")
    print(f"Stacked bar counts CSV: {heatmap.get('stacked_path', '')}")
    print(f"Figures directory: {heatmap.get('figures_dir', '')}")
    print(f"Report: {report}")
    print(f"Manifest: {manifest}")
    failures = [result for result in results if result.get("returncode") != 0]
    if failures:
        print("Failures or missing outputs:")
        for failure in failures:
            print(f"  - {failure.get('parameter_set')}/{failure['method']}: returncode={failure.get('returncode')}")
        return 1 if args.mode == "run" else 0
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

