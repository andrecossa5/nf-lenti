#!/usr/bin/env python3
"""
collapse_outputs.py

Aggregate outputs from nf-lenti bulk pipeline, generate summary statistics, QC plots, and an HTML report.
- Collects per-sample data and metrics.
- Aggregates corrected GBC counts across samples.
- Generates summary CSVs and QC plots.
- Produces a run-level summary JSON and HTML report.
"""

import os
import argparse
import numpy as np
import pandas as pd
import json
from pathlib import Path
import matplotlib.pyplot as plt
import base64
import plotly.graph_objs as go
import math
from matplotlib.ticker import MaxNLocator, ScalarFormatter, NullFormatter

import sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "common"))
from utils import encode_image_b64

def encode_html_b64(path):
    with open(path, "rb") as f:
        return base64.b64encode(f.read()).decode("utf-8")

my_parser = argparse.ArgumentParser(
    prog="collapse_outputs",
    description="Aggregate outputs from nf-lenti bulk pipeline, generate summary statistics, QC plots, and an HTML report."
)

my_parser.add_argument(
    "-i", "--input",
    type=str,
    required=True,
    help="Path to results/perturb_bulk folder."
)

my_parser.add_argument(
    "-o", "--output",
    type=str,
    default="summary",
    help="Name of summary folder to create in the current working directory. Default: 'summary'."
)

my_parser.add_argument(
    "--template",
    type=str,
    required=True,
    help="Path to HTML template."
)

def main():

    args = my_parser.parse_args()
    path_input = Path(args.input).resolve()
    path_output = Path(os.getcwd()) / args.output
    path_output.mkdir(exist_ok=True)

    sample_tables, all_samples = {}, []
    run_info, parameters, metrics_info, pipeline = {}, {}, {}, {}

    # Collect per-sample data from JSON
    for sample_dir in path_input.iterdir():
        if not sample_dir.is_dir() or sample_dir.name == path_output.name:
            continue
        summary_json = sample_dir / "run_summary.json"
        if not summary_json.exists():
            continue
        with open(summary_json) as f:
            stats = json.load(f)
        sample_id = stats.get("sample", sample_dir.name)
        corrected_file = sample_dir / "GBC_counts_corrected.csv"
        raw_file = sample_dir / "GBC_raw_counts.csv.gz"
        if corrected_file.exists() and raw_file.exists():
            corrected = pd.read_csv(corrected_file, index_col=0).rename(columns=lambda c: "corrected_counts")
            raw = pd.read_csv(raw_file, index_col=0).rename(columns=lambda c: "raw_counts")
            df = raw.join(corrected, how="outer").fillna(0).astype(int).assign(sample=sample_id)
            sample_tables[sample_id] = df
        if not run_info: run_info = stats.get("run_info", {})
        if not parameters: parameters = stats.get("parameters", {})
        if not pipeline: pipeline = stats.get("pipeline", {})
        if not metrics_info: metrics_info = stats.get("metrics_info", {})
        qc_images = stats.get("qc_images", [])
        if qc_images and all(isinstance(img, dict) and "title" in img for img in qc_images):
            def find_idx(qc_images, key):
                for i, img in enumerate(qc_images):
                    if key in img["title"].lower():
                        return i
                return -1
            idx_raw = find_idx(qc_images, "raw")
            idx_corr = find_idx(qc_images, "corrected")
            idx_ba = find_idx(qc_images, "before")
            order = [idx_raw, idx_corr, idx_ba]
            if all(i >= 0 for i in order):
                qc_images = [qc_images[i] for i in order]
        all_samples.append({
            "id": sample_id,
            "metrics": stats.get("metrics", {}),
            "qc_images": qc_images
        })

    if not sample_tables:
        raise FileNotFoundError(f"No valid samples found under {path_input}")

    # Aggregate across samples
    corrected_tables = {
        sample_dir.name: pd.read_csv(sample_dir / "GBC_counts_corrected.csv", index_col=0).assign(sample=sample_dir.name)
        for sample_dir in path_input.iterdir()
        if sample_dir.is_dir() and (sample_dir / "GBC_counts_corrected.csv").exists()
    }
    if not corrected_tables:
        raise FileNotFoundError(f"No GBC_counts_corrected.csv files found under {path_input}")

    pd.concat(corrected_tables.values()).to_csv(path_output / "bulk_GBC_reference.csv")

    samples = list(corrected_tables.keys())
    n = len(samples)
    C = np.zeros((n, n))
    for i, x in enumerate(samples):
        for j, y in enumerate(samples):
            if i >= j:
                C[i, j] = C[j, i] = len(set(corrected_tables[x].index) & set(corrected_tables[y].index))
    common_df = pd.DataFrame(C, index=samples, columns=samples)
    common_df.to_csv(path_output / "common.csv")

    # Run-level QC plots (aggregate all samples)
    run_qc_imgs = []
    merged = pd.concat(sample_tables.values())

    # 1. Unique GBC Counts Bar Chart
    samples_order = [s["id"] for s in all_samples]
    raw_gbcs = [s["metrics"].get("n_unique_GBCs", 0) for s in all_samples]
    filtered_gbcs = [s["metrics"].get("n_filtered_GBCs", 0) for s in all_samples]
    corrected_gbcs = [s["metrics"].get("n_corrected_GBCs", 0) for s in all_samples]
    y_all = raw_gbcs + filtered_gbcs + corrected_gbcs
    y_min = min([y for y in y_all if y > 0])
    y_max = max(y_all)
    log_min, log_max = math.floor(np.log10(y_min)), math.ceil(np.log10(y_max))
    n_ticks = min(8, log_max - log_min + 1)
    tick_vals = np.logspace(log_min, log_max, n_ticks) if n_ticks > 1 else [y_min, y_max]
    tick_text = [f"{int(v):,}" for v in tick_vals]
    x = range(len(samples_order))
    width = 0.25
    # Scale figure width with number of samples to prevent label overlap
    fig_width = max(8, len(samples_order) * 0.8)  # Minimum 8 inches, then scale with samples
    plt.figure(figsize=(fig_width, 5))
    plt.bar([i - width for i in x], raw_gbcs, width=width, label="Raw GBCs")
    plt.bar(x, filtered_gbcs, width=width, label="Filtered GBCs")
    plt.bar([i + width for i in x], corrected_gbcs, width=width, label="Corrected GBCs")
    plt.xticks(x, samples_order, rotation=45, ha='right')  # Improved rotation and alignment
    plt.yscale("log")
    plt.ylabel("Unique GBCs")
    plt.title("Unique GBCs Before and After Filtering & Correction")
    plt.legend()
    plt.tight_layout()
    plt.gca().set_yticks(tick_vals)
    plt.gca().set_yticklabels(tick_text)
    plt.gca().yaxis.set_major_formatter(ScalarFormatter())
    plt.gca().yaxis.set_minor_formatter(NullFormatter())
    img_bar = path_output / "run_gbc_summary.png"
    plt.savefig(img_bar, bbox_inches="tight", dpi=300)
    plt.close()

    # Interactive bar chart
    fig_bar = go.Figure([
        go.Bar(x=samples_order, y=raw_gbcs, name="Raw GBCs", marker_color="royalblue"),
        go.Bar(x=samples_order, y=filtered_gbcs, name="Filtered GBCs", marker_color="orange"),
        go.Bar(x=samples_order, y=corrected_gbcs, name="Corrected GBCs", marker_color="green")
    ])
    fig_bar.update_layout(
        barmode="group",
        title=dict(text="Unique GBCs Before and After Filtering & Correction", font=dict(size=16, color="black")),
        xaxis=dict(
            title="Sample", tickangle=15, tickfont=dict(size=16, family="Arial, sans-serif", color="black"),
            showgrid=False, linecolor='black', linewidth=1, mirror=True, showline=True,
            tickvals=samples_order, ticktext=samples_order,
        ),
        yaxis=dict(
            title="Unique GBCs", type="log", tickvals=tick_vals, ticktext=tick_text,
            tickfont=dict(size=16, family="Arial, sans-serif", color="black"),
            showgrid=True, gridcolor="#e5e5e5", linecolor='black', linewidth=1, mirror=True, showline=True,
        ),
        font=dict(family="Arial, sans-serif", size=14, color="black"),
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1, font=dict(size=13, color="black")),
        margin=dict(l=60, r=30, t=60, b=60), plot_bgcolor="white", paper_bgcolor="white"
    )
    interactive_bar_html = path_output / "run_gbc_summary_interactive.html"
    fig_bar.write_html(str(interactive_bar_html), include_plotlyjs="cdn", config={"displayModeBar": False, "displaylogo": False, "responsive": True})
    run_qc_imgs.append({
        "title": "Unique GBCs Before and After Filtering & Correction",
        "image": encode_image_b64(img_bar),
        "interactive_html_b64": encode_html_b64(interactive_bar_html) if interactive_bar_html.exists() else None
    })

    # 2. Cumulative distribution of Corrected GBC Abundance
    corrected_values = merged['corrected_counts'][merged['corrected_counts'] > 0]
    sorted_counts = np.sort(corrected_values)
    cdf = np.arange(1, len(sorted_counts) + 1) / len(sorted_counts)
    plt.figure(figsize=(6, 4))
    plt.plot(sorted_counts, cdf, marker='.', linestyle='none', color="green", alpha=0.7)
    plt.xscale('log')
    plt.xlabel("Corrected GBC Read Count (log10 scale)")
    plt.ylabel("Cumulative Fraction of GBCs")
    plt.title("Cumulative Distribution of Corrected GBC Abundance")
    plt.grid(True, which="both", ls="--", lw=0.5, alpha=0.6)
    plt.tight_layout()
    img_cdf = path_output / "run_corrected_abundance_cdf.png"
    plt.savefig(img_cdf, bbox_inches="tight", dpi=300)
    plt.close()
    fig_cum = go.Figure([
        go.Scatter(x=sorted_counts, y=cdf, mode='markers', marker=dict(color="green", size=5, opacity=0.7), name="CDF")
    ])
    fig_cum.update_layout(
        title=dict(text="Cumulative Distribution of Corrected GBC Abundance", font=dict(size=16, color="black")),
        xaxis=dict(
            title="Corrected GBC Read Count (log10 scale)", type="log",
            tickfont=dict(size=16, family="Arial, sans-serif", color="black"),
            showgrid=True, gridcolor="#e5e5e5", linecolor='black', linewidth=1, mirror=True, showline=True,
        ),
        yaxis=dict(
            title="Cumulative Fraction of GBCs",
            tickfont=dict(size=16, family="Arial, sans-serif", color="black"),
            showgrid=True, gridcolor="#e5e5e5", linecolor='black', linewidth=1, mirror=True, showline=True,
            range=[0, 1]
        ),
        font=dict(family="Arial, sans-serif", size=14, color="black"),
        margin=dict(l=60, r=30, t=60, b=60),
        plot_bgcolor="white", paper_bgcolor="white"
    )
    interactive_cdf_html = path_output / "run_corrected_abundance_cdf_interactive.html"
    fig_cum.write_html(str(interactive_cdf_html), include_plotlyjs="cdn", config={"displayModeBar": False, "displaylogo": False, "responsive": True})
    run_qc_imgs.append({
        "title": "Cumulative Distribution of Corrected GBC Abundance",
        "image": encode_image_b64(img_cdf),
        "interactive_html_b64": encode_html_b64(interactive_cdf_html) if interactive_cdf_html.exists() else None
    })

    # 3. Prevalence of Corrected GBCs Across Samples
    corrected_matrix = pd.DataFrame({
        sample: tbl['corrected_counts'] if 'corrected_counts' in tbl else pd.Series(dtype=int)
        for sample, tbl in sample_tables.items()
    }).fillna(0).astype(int)
    prevalence = (corrected_matrix > 0).sum(axis=1)
    plt.figure(figsize=(6,4))
    plt.hist(prevalence, bins=np.arange(1, len(samples)+2)-0.5, color="purple", alpha=0.7)
    plt.xlabel("Number of Samples GBC is Present In")
    plt.ylabel("Number of GBCs")
    plt.title("Prevalence of Corrected GBCs Across Samples")
    plt.xticks(range(1, len(samples)+1))
    plt.gca().yaxis.set_major_locator(MaxNLocator(nbins=8, integer=True))
    plt.tight_layout()
    img_prev = path_output / "run_gbc_prevalence.png"
    plt.savefig(img_prev, bbox_inches="tight", dpi=300)
    plt.close()
    interactive_prev_html = path_output / "run_gbc_prevalence_interactive.html"
    bin_edges = np.arange(1, len(samples)+2) - 0.5
    hist, bins = np.histogram(prevalence, bins=bin_edges)
    bin_centers = (bins[:-1] + bins[1:]) / 2
    fig_prev = go.Figure([
        go.Bar(
            x=bin_centers, y=hist, marker_color='purple',
            hovertemplate='Samples: %{x}<br>GBCs: %{y}<extra></extra>', name="Prevalence"
        )
    ])
    fig_prev.update_layout(
        title=dict(text="Prevalence of Corrected GBCs Across Samples", font=dict(size=16, color="black")),
        xaxis=dict(
            title="Number of Samples GBC is Present In",
            tickmode="array", tickvals=list(range(1, len(samples)+1)),
            ticktext=[str(i) for i in range(1, len(samples)+1)],
            tickfont=dict(size=16, family="Arial, sans-serif", color="black"),
            showgrid=False, zeroline=False, linecolor='black', linewidth=1, mirror=True,
        ),
        yaxis=dict(
            title="Number of GBCs",
            tickfont=dict(size=16, family="Arial, sans-serif", color="black"),
            showgrid=True, gridcolor="#e5e5e5", zeroline=False, linecolor='black', linewidth=1, mirror=True,
        ),
        font=dict(family="Arial, sans-serif", size=14, color="black"),
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1, font=dict(size=13, color="black")),
        margin=dict(l=60, r=30, t=60, b=60),
        plot_bgcolor="white", paper_bgcolor="white"
    )
    fig_prev.update_xaxes(
        tickfont=dict(size=16, family="Arial, sans-serif", color="black"),
        title=dict(font=dict(size=18, color="black"))
    )
    fig_prev.update_yaxes(
        tickfont=dict(size=16, family="Arial, sans-serif", color="black"),
        title=dict(font=dict(size=18, color="black"))
    )
    fig_prev.write_html(str(interactive_prev_html), include_plotlyjs="cdn", config={"displayModeBar": False, "displaylogo": False, "responsive": True})
    run_qc_imgs.append({
        "title": "Prevalence of Corrected GBCs Across Samples",
        "image": encode_image_b64(img_prev),
        "interactive_html_b64": encode_html_b64(interactive_prev_html) if interactive_prev_html.exists() else None
    })

    # 4. Heatmap of common clones between samples
    fig, ax = plt.subplots(figsize=(max(6, n), max(6, n)))
    cmap = plt.get_cmap("mako") if "mako" in plt.colormaps() else plt.get_cmap("viridis")
    im = ax.imshow(common_df.values, cmap=cmap, vmin=0, vmax=np.percentile(common_df.values, 99))
    ax.set_xticks(np.arange(n))
    ax.set_yticks(np.arange(n))
    ax.set_xticklabels(samples, rotation=90, fontsize=8)
    ax.set_yticklabels(samples, fontsize=8)
    ax.set_title("Number of Common clones")
    norm = im.norm
    for i in range(n):
        for j in range(n):
            value = common_df.values[i, j]
            normed = norm(value)
            r, g, b, _ = cmap(normed)
            brightness = (r*299 + g*587 + b*114) / 1000
            text_color = "black" if brightness > 0.5 else "white"
            ax.text(j, i, f"{int(value)}", ha="center", va="center", color=text_color, fontsize=10)
    fig.colorbar(im, ax=ax, fraction=0.05, aspect=35, pad=0.02, label="n common clones")
    fig.tight_layout()
    img_heatmap = path_output / "common_heatmap.png"
    fig.savefig(img_heatmap, dpi=500)
    plt.close(fig)
    run_qc_imgs.append({
        "title": "Common clones Heatmap",
        "image": encode_image_b64(img_heatmap)
    })

    # Final run_summary.json
    run_summary = {
        "run_type": "bulk",
        "pipeline": pipeline,
        "run_info": run_info,
        "date": str(pd.Timestamp.today().date()),
        "parameters": parameters,
        "metrics_info": metrics_info,
        "run_level_qc": run_qc_imgs,
        "samples": all_samples
    }
    with open(path_output / "run_summary.json", "w") as f:
        json.dump(run_summary, f, indent=2)

    # Generate HTML
    with open(args.template) as f:
        template = f.read()
    html_out = template.replace("{{REPORT_JSON}}", json.dumps(run_summary))
    with open(path_output / "run_report.html", "w") as f:
        f.write(html_out)
    print(f"[collapse_outputs] Summary + HTML written to: {path_output}")

if __name__ == "__main__":
    main()