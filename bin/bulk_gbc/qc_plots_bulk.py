#!/usr/bin/env python3

"""
qc_plots_bulk.py

Generate QC plots for Bulk GBC workflow (no seaborn).
- Plots the distribution of raw and corrected GBC read counts.
- Plots before vs after correction for GBC read counts.
- Optionally plots spike-in control trends if provided.
"""

import os
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

# ----------------------------------------------------
# Argument parsing
# ----------------------------------------------------
parser = argparse.ArgumentParser(
    prog="qc_plots_bulk",
    description="Generate QC plots for Bulk GBC workflow (no seaborn)."
)

parser.add_argument(
    "--raw_counts",
    type=str,
    required=True,
    help="Path to raw counts CSV."
)
parser.add_argument(
    "--corrected_counts",
    type=str,
    required=True,
    help="Path to corrected counts CSV."
)
parser.add_argument(
    "--correction_df",
    type=str,
    required=True,
    help="Path to correction dataframe CSV."
)
parser.add_argument(
    "--outdir",
    type=str,
    required=True,
    help="Output directory for plots."
)
parser.add_argument(
    "--spikeins",
    type=str,
    default=None,
    help="Path to spike-ins CSV (optional)."
)
parser.add_argument(
    "--n_reads",
    type=int,
    default=1000,
    help="Minimum number of reads (default: 1000)."
)

args = parser.parse_args()

os.makedirs(args.outdir, exist_ok=True)

# ----------------------------------------------------
# Helpers
# ----------------------------------------------------
def kdeplot(data, ax, label=None, color=None):
    """Simple KDE with matplotlib."""
    kde = gaussian_kde(data)
    x = np.linspace(min(data), max(data), 200)
    y = kde(x)
    ax.plot(x, y, label=label, color=color)
    ax.fill_between(x, 0, y, alpha=0.3, color=color)

def scatter_with_fit(x, y, ax, xlabel, ylabel, title):
    """Scatter + simple linear fit."""
    ax.scatter(x, y, color="black", s=10)
    if len(x) > 1:
        m, b = np.polyfit(x, y, 1)
        ax.plot(x, m*x + b, color="red", lw=1)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

def standardize_read_count(df):
    for col in ["read_counts", "count_reads", "count_read"]:
        if col in df.columns:
            df = df.rename(columns={col: "read_count"})
    return df

# ----------------------------------------------------
# Load data
# ----------------------------------------------------
raw_counts = standardize_read_count(pd.read_csv(args.raw_counts, index_col=0))
corrected_counts = standardize_read_count(pd.read_csv(args.corrected_counts, index_col=0))
correction_df = pd.read_csv(args.correction_df)

# ----------------------------------------------------
# Plot 1: Raw counts distribution
# ----------------------------------------------------
fig, ax = plt.subplots(figsize=(6,4))
log_counts = np.log10(raw_counts["read_count"] + 1)
kdeplot(log_counts, ax, label="Raw counts", color="blue")
ax.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.7)
median_val = np.median(log_counts)
ax.axvline(median_val, color='red', linestyle='--', linewidth=1.2, label=f"Median: {median_val:.2f}")
ax.set_title("Distribution of Raw GBC Read Counts (log10 scale)")
ax.set_xlabel("log10(Read Count + 1)")
ax.set_ylabel("Density")
ax.legend()
fig.tight_layout()
fig.savefig(os.path.join(args.outdir, "raw_counts_distribution.png"), dpi=300)
plt.close(fig)

# ----------------------------------------------------
# Plot 2: Corrected counts distribution
# ----------------------------------------------------
fig, ax = plt.subplots(figsize=(6,4))
log_counts_corr = np.log10(corrected_counts["read_count"] + 1)
kdeplot(log_counts_corr, ax, label="Corrected counts", color="green")
ax.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.7)
median_val_corr = np.median(log_counts_corr)
ax.axvline(median_val_corr, color='red', linestyle='--', linewidth=1.2, label=f"Median: {median_val_corr:.2f}")
ax.set_title("Distribution of Corrected GBC Read Counts (log10 scale)")
ax.set_xlabel("log10(Read Count + 1)")
ax.set_ylabel("Density")
ax.legend()
fig.tight_layout()
fig.savefig(os.path.join(args.outdir, "corrected_counts_distribution.png"), dpi=300)
plt.close(fig)

# ----------------------------------------------------
# Plot 3: Reads before vs after correction
# ----------------------------------------------------
fig, ax = plt.subplots(figsize=(5,5))
before = raw_counts.loc[corrected_counts.index, "read_count"]
after = corrected_counts["read_count"]
scatter_with_fit(
    np.log10(before+1), 
    np.log10(after+1), 
    ax,
    xlabel="log10 Reads Before Correction",
    ylabel="log10 Reads After Correction",
    title="GBC Read Counts: Before vs After Correction"
)
fig.tight_layout()
fig.savefig(os.path.join(args.outdir, "before_vs_after_correction.png"), dpi=300)
plt.close(fig)

# ----------------------------------------------------
# Spike-ins plot (if provided)
# ----------------------------------------------------
if args.spikeins and os.path.exists(args.spikeins):
    spikes = pd.read_csv(args.spikeins)
    numeric_cols = spikes.select_dtypes(include=[np.number]).columns

    if len(numeric_cols) >= 2:
        xcol, ycol = numeric_cols[:2]
        fig, ax = plt.subplots(figsize=(5,5))
        scatter_with_fit(
            np.log10(spikes[xcol] + 1),
            np.log10(spikes[ycol] + 1),
            ax,
            xlabel=f"log10 {xcol} + 1",
            ylabel=f"log10 {ycol} + 1",
            title="Spike-in Control Trend"
        )
        fig.tight_layout()
        fig.savefig(os.path.join(args.outdir, "Spike-in Control Trend.png"), dpi=300)
        plt.close(fig)
        print(f"[qc_plots_bulk] Spike-ins plot using {xcol} vs {ycol}")

else:
    print("[qc_plots_bulk] No spike-ins provided.")
