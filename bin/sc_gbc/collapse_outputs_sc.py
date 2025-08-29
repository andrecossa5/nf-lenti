#!/usr/bin/env python3
import os
import json
import argparse
from pathlib import Path
from datetime import date
import base64

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.graph_objs as go
import numpy as np
import math
import matplotlib.ticker as mticker

import sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "common"))
from utils import encode_plot_to_base64, get_sample_id_from_json, get_metrics_for_sample

"""
collapse_outputs_sc.py

Collapse single-cell GBC outputs into a run-level summary and HTML report.
"""

def get_figsize(n_samples, base_width=12, per_sample=0.7, height=4):
    return (max(base_width, n_samples * per_sample), height)

def nice_step(x):
    for s in [0.1, 0.2, 0.25, 0.5, 1, 2, 5, 10]:
        if x <= s:
            return s
    return 1

def log_ticks(y_min, y_max, n_ticks=8):
    if y_min > 0 and y_max > 0 and y_max > y_min:
        log_min = math.floor(np.log10(y_min))
        log_max = math.ceil(np.log10(y_max))
        ticks = [y_min]
        for exp in range(log_min, log_max + 1):
            val = 10 ** exp
            if y_min < val < y_max:
                ticks.append(val)
                mid = 5 * 10 ** exp
                if y_min < mid < y_max:
                    ticks.append(mid)
        if y_max not in ticks:
            ticks.append(y_max)
        return sorted(set(int(round(y)) for y in ticks))
    else:
        return [y_min, y_max] if y_min != y_max else [y_min]

def main():
    parser = argparse.ArgumentParser(
        prog="collapse_outputs_sc",
        description="Collapse single-cell GBC outputs into a run-level summary and HTML report."
    )
    parser.add_argument("-i", "--input", type=str, required=True,
                        help="Path to input sc_outdir (parent directory containing per-sample subdirectories).")
    parser.add_argument("--template", type=str, required=True,
                        help="Path to HTML template for the report.")
    parser.add_argument("-o", "--output", type=str, default="sc_summary",
                        help="Output directory for summary files. Default: 'sc_summary'.")
    args = parser.parse_args()

    path_input = Path(args.input).resolve()
    path_output = Path(os.getcwd()) / args.output
    path_output.mkdir(exist_ok=True)

    samples_collected, umi_data, moi_data, reads_data = [], [], [], []
    run_info = parameters = pipeline = metrics_info = None
    matched_dirs = 0

    print(f"[INFO] Scanning {path_input} for per-sample outputs...")

    for root, _, files in os.walk(path_input):
        root = Path(root)
        if not ("run_summary.json" in files and "CBC_GBC_combos.tsv.gz" in files):
            continue
        matched_dirs += 1
        json_path = root / "run_summary.json"
        combos_path = root / "CBC_GBC_combos.tsv.gz"
        try:
            with open(json_path) as f:
                data = json.load(f)
        except Exception as e:
            print(f"[WARNING] Cannot read run_summary.json in {root}: {e}")
            continue
        if run_info is None: run_info = data.get("run_info", {})
        if parameters is None: parameters = data.get("parameters", {})
        if pipeline is None: pipeline = data.get("pipeline", {})
        if metrics_info is None: metrics_info = data.get("metrics_info", {})
        sample_id = get_sample_id_from_json(data, fallback=root.name)
        m = get_metrics_for_sample(data, sample_id)
        if "samples" in data:
            for s in data["samples"]:
                if not s.get("id"):
                    s["id"] = sample_id
            samples_collected.extend(data["samples"])
        reads = m.get("total_reads") or m.get("mapped_reads")
        if reads is not None:
            reads_data.append({"sample": sample_id, "reads": reads})
        moi = m.get("median_moi") or m.get("observed_moi_median")
        if moi is not None:
            moi_data.append({"sample": sample_id, "moi": moi})
        try:
            df_umi = pd.read_csv(combos_path, sep="\t", compression="gzip")
        except Exception as e:
            print(f"[WARNING] Failed to read UMIs from {combos_path}: {e}")
            continue
        if df_umi is None or df_umi.empty or "umi" not in df_umi.columns:
            continue
        df_umi["umi"] = pd.to_numeric(df_umi["umi"], errors="coerce")
        df_umi = df_umi.dropna(subset=["umi"])
        df_umi = df_umi[df_umi["umi"] > 0]
        if df_umi.empty:
            continue
        umi_data.extend([{"sample": sample_id, "UMIs": int(u)} for u in df_umi["umi"]])

    if matched_dirs == 0:
        raise RuntimeError("No valid sample directories found.")
    if not samples_collected:
        raise RuntimeError("Found directories but no samples collected.")

    run_level_qc = []

    # Total Reads
    if reads_data:
        df_reads = pd.DataFrame(reads_data)
        fig, ax = plt.subplots(figsize=get_figsize(len(df_reads)))
        sns.barplot(x="sample", y="reads", data=df_reads, palette="mako", ax=ax)
        ax.set_title("Total Reads per Sample", fontsize=20, color="black")
        ax.set_xlabel("Sample", fontsize=17, color="black")
        ax.set_ylabel("Number of Reads", fontsize=17, color="black")
        ax.set_yscale("log")
        y_min, y_max = df_reads["reads"].min(), df_reads["reads"].max()
        log_min, log_max = math.floor(np.log10(y_min)), math.ceil(np.log10(y_max))
        n_ticks = min(8, log_max - log_min + 1)
        yticks = np.logspace(log_min, log_max, n_ticks) if n_ticks > 1 else [y_min, y_max]
        ax.set_yticks(yticks)
        ax.yaxis.set_major_formatter(mticker.ScalarFormatter())
        ax.yaxis.set_minor_formatter(mticker.NullFormatter())
        ax.set_xticklabels(df_reads["sample"], rotation=45, ha='right', fontsize=15, color="black")
        ax.tick_params(axis='y', labelsize=15, colors="black")
        ax.tick_params(axis='x', labelsize=15, colors="black")
        fig.tight_layout()
        fig.savefig(path_output / "total_reads_per_sample.png", dpi=300)
        fig_reads = go.Figure([go.Bar(
            x=df_reads["sample"], y=df_reads["reads"], marker_color="teal"
        )])
        tick_text = [f"{int(y):,}" for y in yticks]
        fig_reads.update_layout(
            title=dict(text="Total Reads per Sample", font=dict(size=24, color="black")),
            xaxis=dict(
                title="Sample", ticks="outside", tickangle=15, tickfont=dict(size=18, color="black"),
                linecolor='black', linewidth=1, mirror=True, showgrid=True,
                tickvals=df_reads["sample"].tolist(), ticktext=df_reads["sample"].tolist(), showline=True,
            ),
            yaxis=dict(
                title="Number of Reads", type="log", ticks="outside",
                tickvals=list(yticks), ticktext=tick_text,
                tickfont=dict(size=18, color="black"), linecolor='black', linewidth=1,
                mirror=True, showgrid=True, showline=True,
            ),
            font=dict(family="Arial, sans-serif", size=16, color="black"),
            margin=dict(l=60, r=30, t=60, b=60),
            plot_bgcolor="white", paper_bgcolor="white", height=420, autosize=True,
        )
        interactive_bar_html = path_output / "total_reads_per_sample_interactive.html"
        fig_reads.write_html(
            str(interactive_bar_html),
            include_plotlyjs="cdn",
            config={"displayModeBar": False, "displaylogo": False, "responsive": True}
        )
        interactive_html_b64 = None
        if interactive_bar_html.exists():
            with open(interactive_bar_html, "rb") as f:
                interactive_html_b64 = base64.b64encode(f.read()).decode("utf-8")
        run_level_qc.append({
            "title": "Total Reads per Sample",
            "image": encode_plot_to_base64(fig),
            "interactive_html_b64": interactive_html_b64
        })
        plt.close(fig)

    # Median MOI
    if moi_data:
        df_moi = pd.DataFrame(moi_data)
        fig, ax = plt.subplots(figsize=get_figsize(len(df_moi)))
        sns.barplot(x="sample", y="moi", data=df_moi, palette="crest", ax=ax)
        ax.set_title("Median MOI per Sample", fontsize=20, color="black")
        ax.set_xlabel("Sample", fontsize=17, color="black")
        ax.set_ylabel("Median MOI", fontsize=17, color="black")
        moi_min, moi_max = float(df_moi["moi"].min()), float(df_moi["moi"].max())
        n_ticks = 8
        if moi_max > moi_min:
            step = nice_step((moi_max - moi_min) / (n_ticks - 1))
            start = np.floor(moi_min / step) * step
            end = np.ceil(moi_max / step) * step
            yticks = np.round(np.arange(start, end + step/2, step), 6)
        else:
            yticks = np.array([moi_min])
        ax.set_yticks(yticks)
        ax.yaxis.set_major_formatter(mticker.FormatStrFormatter('%.2f'))
        ax.yaxis.set_minor_formatter(mticker.NullFormatter())
        ax.set_xticklabels(df_moi["sample"], rotation=45, ha='right', fontsize=15, color="black")
        ax.tick_params(axis='y', labelsize=15, colors="black", direction='out', length=6, width=1)
        ax.tick_params(axis='x', labelsize=15, colors="black", direction='out', length=6, width=1)
        fig.tight_layout()
        fig.savefig(path_output / "median_moi_per_sample.png", dpi=300)
        moi_vals = df_moi["moi"].values
        y_min, y_max = float(np.min(moi_vals)), float(np.max(moi_vals))
        if y_max > y_min:
            step = nice_step((y_max - y_min) / (n_ticks - 1))
            start = np.floor(y_min / step) * step
            end = np.ceil(y_max / step) * step
            yticks_interactive = np.round(np.arange(start, end + step/2, step), 6)
            yticktext = [f"{y:.2f}" for y in yticks_interactive]
        else:
            yticks_interactive = [y_min]
            yticktext = [f"{y_min:.2f}"]
        fig_moi = go.Figure([go.Bar(x=df_moi["sample"], y=df_moi["moi"], marker_color="darkcyan")])
        fig_moi.update_layout(
            title=dict(text="Median MOI per Sample", font=dict(size=24, color="black")),
            xaxis=dict(
                title="Sample", ticks="outside", tickangle=45, tickfont=dict(size=18, color="black"),
                linecolor='black', linewidth=1, mirror=True, showgrid=True,
                tickvals=df_moi["sample"].tolist(), ticktext=df_moi["sample"].tolist(), showline=True,
            ),
            yaxis=dict(
                title="Median MOI", ticks="outside",
                tickvals=list(yticks_interactive), ticktext=yticktext,
                tickfont=dict(size=18, color="black"), linecolor='black', linewidth=1,
                mirror=True, showgrid=True, showline=True,
            ),
            font=dict(family="Arial, sans-serif", size=16, color="black"),
            margin=dict(l=60, r=30, t=60, b=60),
            plot_bgcolor="white", paper_bgcolor="white", height=420, autosize=True
        )
        interactive_moi_html = path_output / "median_moi_per_sample_interactive.html"
        fig_moi.write_html(
            str(interactive_moi_html),
            include_plotlyjs="cdn",
            config={"displayModeBar": False, "displaylogo": False, "responsive": True}
        )
        interactive_html_b64 = None
        if interactive_moi_html.exists():
            with open(interactive_moi_html, "rb") as f:
                interactive_html_b64 = base64.b64encode(f.read()).decode("utf-8")
        run_level_qc.append({
            "title": "Median MOI per Sample",
            "image": encode_plot_to_base64(fig),
            "interactive_html_b64": interactive_html_b64
        })
        plt.close(fig)

    # UMI Distribution
    if umi_data:
        df_umi_plot = pd.DataFrame(umi_data)
        fig, ax = plt.subplots(figsize=get_figsize(len(df_umi_plot["sample"].unique())))
        sns.boxplot(x="sample", y="UMIs", data=df_umi_plot, palette="viridis", ax=ax)
        ax.set_title("UMI Distribution per Sample", fontsize=20, color="black")
        ax.set_xlabel("Sample", fontsize=17, color="black")
        ax.set_ylabel("Number of UMIs per CBC-GBC", fontsize=17, color="black")
        ax.set_yscale("log")
        y_min, y_max = df_umi_plot["UMIs"].min(), df_umi_plot["UMIs"].max()
        yticks = log_ticks(y_min, y_max)
        ax.set_yticks(yticks)
        ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: f"{int(x):,}" if x >= 1 else f"{x:.2f}"))
        ax.yaxis.set_minor_formatter(mticker.NullFormatter())
        ax.set_xticklabels(df_umi_plot["sample"].unique(), rotation=45, ha='right', fontsize=15, color="black")
        ax.tick_params(axis='y', labelsize=15, colors="black")
        ax.tick_params(axis='x', labelsize=15, colors="black")
        fig.tight_layout()
        fig.savefig(path_output / "umi_distribution_per_sample.png", dpi=300)
        fig_umi = go.Figure()
        unique_samples = list(df_umi_plot["sample"].unique())
        for sample in unique_samples:
            sample_umis = df_umi_plot[df_umi_plot["sample"] == sample]["UMIs"]
            fig_umi.add_trace(go.Box(y=sample_umis, name=sample, boxpoints='outliers', marker_color="seagreen"))
        tick_text = [f"{int(y):,}" if y >= 1 else f"{y:.2f}" for y in yticks]
        fig_umi.update_layout(
            title=dict(text="UMI Distribution per Sample", font=dict(size=24, color="black")),
            xaxis=dict(
                title="Sample", tickangle=45, tickfont=dict(size=18, color="black"),
                linecolor='black', linewidth=1, mirror=True, showgrid=True,
                tickvals=unique_samples, ticktext=unique_samples, ticks="outside", showline=True,
            ),
            yaxis=dict(
                title="Number of UMIs per CBC-GBC", type="log",
                tickvals=yticks, ticktext=tick_text,
                tickfont=dict(size=18, color="black"), linecolor='black', linewidth=1,
                mirror=True, showgrid=True, ticks="outside", showline=True,
            ),
            font=dict(family="Arial, sans-serif", size=16, color="black"),
            margin=dict(l=60, r=30, t=60, b=60),
            plot_bgcolor="white", paper_bgcolor="white", height=420, autosize=True,
        )
        interactive_umi_html = path_output / "umi_distribution_per_sample_interactive.html"
        fig_umi.write_html(
            str(interactive_umi_html),
            include_plotlyjs="cdn",
            config={"displayModeBar": False, "displaylogo": False, "responsive": True}
        )
        interactive_html_b64 = None
        if interactive_umi_html.exists():
            with open(interactive_umi_html, "rb") as f:
                interactive_html_b64 = base64.b64encode(f.read()).decode("utf-8")
        run_level_qc.append({
            "title": "UMI Distribution per Sample",
            "image": encode_plot_to_base64(fig),
            "interactive_html_b64": interactive_html_b64
        })
        plt.close(fig)

    sample_ids = [s.get("id") for s in samples_collected if s.get("id")]
    run_summary = {
        "id": ",".join(sample_ids),
        "run_type": "sc",
        "date": date.today().isoformat(),
        "run_info": run_info if run_info else {},
        "pipeline": pipeline if pipeline else {},
        "metrics_info": metrics_info if metrics_info else {},
        "parameters": parameters if parameters else {},
        "run_level_qc": run_level_qc,
        "samples": samples_collected
    }

    with open(path_output / "run_summary.json", "w") as f:
        json.dump(run_summary, f, indent=2)

    with open(args.template) as f:
        template = f.read()
    html_out = template.replace("{{REPORT_JSON}}", json.dumps(run_summary))
    with open(path_output / "run_report.html", "w") as f:
        f.write(html_out)

    print(f"[INFO] Wrote {path_output/'run_summary.json'} and run_report.html")

if __name__ == "__main__":
    main()