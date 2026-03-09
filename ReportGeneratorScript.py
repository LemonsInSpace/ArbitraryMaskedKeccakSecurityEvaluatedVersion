#!/usr/bin/env python3

import sys
import json
import base64
import pathlib
import numpy as np
import matplotlib.pyplot as plt
from io import BytesIO

# ============================================================
# CONFIG
# ============================================================

TVLA_THRESHOLD = 4.5
ORDERS = ["t1", "t_var", "t2", "t2_group", "t3", "t3_group", "t4", "t4_group"]

# ============================================================
# Utilities
# ============================================================

def load_npy(path):
    if path.exists():
        return np.load(path)
    return None

def max_abs(arr):
    if arr is None:
        return None, None
    idx = int(np.nanargmax(np.abs(arr)))
    return float(np.abs(arr[idx])), idx

def count_exceed(arr, threshold):
    if arr is None:
        return 0
    return int(np.sum(np.abs(arr) > threshold))

def plot_array(arr, title):
    if arr is None:
        return None

    fig = plt.figure(figsize=(12, 4))
    plt.plot(arr)
    plt.axhline(TVLA_THRESHOLD, linestyle="--")
    plt.axhline(-TVLA_THRESHOLD, linestyle="--")
    plt.title(title)
    plt.tight_layout()

    buf = BytesIO()
    plt.savefig(buf, format="png")
    plt.close(fig)
    buf.seek(0)

    return base64.b64encode(buf.read()).decode("ascii")

# ============================================================
# Main
# ============================================================

if len(sys.argv) != 2:
    print("Usage: python3 view_capture_report.py <capture_directory>")
    sys.exit(1)

capture_dir = pathlib.Path(sys.argv[1])
if not capture_dir.exists():
    print("Directory does not exist")
    sys.exit(1)

print(f"Generating report for: {capture_dir}")

report = {}
plots = {}

# ------------------------------------------------------------
# Load metadata
# ------------------------------------------------------------

meta = None
meta_path = capture_dir / "meta.json"
if meta_path.exists():
    meta = json.loads(meta_path.read_text())
report["meta"] = meta

# ------------------------------------------------------------
# Load TVLA arrays
# ------------------------------------------------------------

tvla_summary = {}

for name in ORDERS:
    arr = load_npy(capture_dir / f"{name}.npy")
    peak, idx = max_abs(arr)
    exceed = count_exceed(arr, TVLA_THRESHOLD)

    tvla_summary[name] = {
        "max_abs": peak,
        "peak_index": idx,
        "exceed_count": exceed,
    }

    plots[name] = plot_array(arr, f"{name} (TVLA)")

report["tvla"] = tvla_summary

# ------------------------------------------------------------
# Mean + variance stats
# ------------------------------------------------------------

delta_mu = load_npy(capture_dir / "delta_mu.npy")
var_fixed = load_npy(capture_dir / "var_fixed.npy")
var_random = load_npy(capture_dir / "var_random.npy")
pooled_var = load_npy(capture_dir / "pooled_var.npy")

report["means"] = {
    "delta_mu_max_abs": max_abs(delta_mu)[0] if delta_mu is not None else None
}

report["variance"] = {
    "var_fixed_mean": float(np.mean(var_fixed)) if var_fixed is not None else None,
    "var_random_mean": float(np.mean(var_random)) if var_random is not None else None,
    "pooled_var_mean": float(np.mean(pooled_var)) if pooled_var is not None else None,
}

plots["delta_mu"] = plot_array(delta_mu, "Delta Mean")
plots["var_fixed"] = plot_array(var_fixed, "Variance Fixed")
plots["var_random"] = plot_array(var_random, "Variance Random")

# ------------------------------------------------------------
# Drift (global mean/variance over traces)
# ------------------------------------------------------------

global_mean_series = load_npy(capture_dir / "global_mean_series.npy")
global_var_series  = load_npy(capture_dir / "global_var_series.npy")

plots["global_mean_series"] = plot_array(global_mean_series, "Global Mean Drift")
plots["global_var_series"]  = plot_array(global_var_series, "Global Variance Drift")

# ------------------------------------------------------------
# Energy stats
# ------------------------------------------------------------

energy_stats_path = capture_dir / "trace_energy_stats.json"
if energy_stats_path.exists():
    report["energy"] = json.loads(energy_stats_path.read_text())

# ------------------------------------------------------------
# Timing stats
# ------------------------------------------------------------

timing_path = capture_dir / "timing_stats.json"
if timing_path.exists():
    report["timing"] = json.loads(timing_path.read_text())

# ------------------------------------------------------------
# Write JSON summary
# ------------------------------------------------------------

json_out = capture_dir / "report.json"
json_out.write_text(json.dumps(report, indent=2))
print(f"Saved: {json_out}")

# ------------------------------------------------------------
# Write HTML report (self-contained)
# ------------------------------------------------------------

html_parts = []
html_parts.append("<html><head><title>TVLA Report</title></head><body>")
html_parts.append(f"<h1>Capture Report: {capture_dir.name}</h1>")

html_parts.append("<h2>Meta</h2>")
html_parts.append(f"<pre>{json.dumps(meta, indent=2)}</pre>")

html_parts.append("<h2>TVLA Summary</h2>")
html_parts.append(f"<pre>{json.dumps(tvla_summary, indent=2)}</pre>")

for name, img in plots.items():
    if img:
        html_parts.append(f"<h3>{name}</h3>")
        html_parts.append(f'<img src="data:image/png;base64,{img}"/>')

html_parts.append("</body></html>")

html_out = capture_dir / "report.html"
html_out.write_text("\n".join(html_parts))

print(f"Saved: {html_out}")
print("Done.")
