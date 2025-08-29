#!/usr/bin/env python3

"""
create_run_summary.py

Create run summary for bulk GBC analysis.
- Loads raw, corrected, and correction dataframe CSVs for a sample.
- Computes key metrics about GBC filtering and correction.
- Encodes QC images for embedding in reports.
- Outputs a summary JSON with pipeline, run info, parameters, metrics, and QC images.
"""

import os
import sys
import argparse
import json
import datetime
import numpy as np
import pandas as pd
from pathlib import Path

sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "common"))
from utils import encode_image_b64, resolve_working_dir_and_user

def main():
    parser = argparse.ArgumentParser(
        prog="create_run_summary",
        description="Create run summary for bulk GBC analysis."
    )
    parser.add_argument(
        "--indir",
        type=str,
        required=True,
        help="Input directory."
    )
    parser.add_argument(
        "--outdir",
        type=str,
        required=True,
        help="Output directory."
    )
    parser.add_argument(
        "--params_outdir",
        type=str,
        required=True,
        default=None,
        help="Original pipeline outdir parameter (for reporting only)."
    )
    parser.add_argument(
        "--anchor_sequence",
        type=str,
        required=True,
        help="Anchor sequence used for GBC extraction."
    )
    parser.add_argument(
        "--sample",
        type=str,
        required=True,
        help="Sample name."
    )
    parser.add_argument(
        "--raw_counts",
        type=str,
        required=True,
        help="Path to raw counts CSV file (before filtering/correction)."
    )
    parser.add_argument(
        "--corrected_counts",
        type=str,
        required=True,
        help="Path to corrected counts CSV file (after filtering/correction)."
    )
    parser.add_argument(
        "--correction_df",
        type=str,
        required=True,
        help="Path to correction dataframe CSV file (mapping degenerate to corrected GBCs)."
    )
    parser.add_argument(
        "--min_n_reads",
        type=int,
        required=False,
        default=1000,
        help="Minimum number of reads that a 'correct' GBC must have to be considered for entry into the final whitelist. Default: 1000."
    )
    parser.add_argument(
        "--hamming_treshold",
        type=int,
        required=False,
        default=3,
        help="Maximum Hamming distance for which sequences can be clustered together."
    )
    parser.add_argument(
        "--qc_images",
        type=str,
        nargs='+',
        required=True,
        help="List of QC image paths."
    )
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    working_dir, user_name = resolve_working_dir_and_user()

    # Load data
    df_raw = pd.read_csv(args.raw_counts, index_col=0)
    df_corrected = pd.read_csv(args.corrected_counts, index_col=0)
    df_correction = pd.read_csv(args.correction_df)
    if 'read_counts' in df_raw.columns:
        df_raw.rename(columns={'read_counts': 'read_count'}, inplace=True)

    # Compute stats
    total_GBC_reads = df_raw['read_count'].sum()
    n_unique_GBCs = df_raw.shape[0]
    n_filtered_GBCs = df_corrected.shape[0]
    n_corrected = df_correction['correct'].nunique() if 'correct' in df_correction.columns else 0
    n_degenerated = df_correction['degenerated'].nunique() if 'degenerated' in df_correction.columns else 0
    median_n_reads_added = float(df_correction['n_reads_degenerated'].median()) if 'n_reads_degenerated' in df_correction.columns else 0.0

    # Encode QC images with titles
    qc_encoded = [
        {"title": Path(p).stem, "image": encode_image_b64(p)}
        for p in args.qc_images if Path(p).exists()
    ]

    # Build summary JSON
    summary_data = {
        "pipeline": {
            "name": "nf-lenti",
            "version": "1.0.0",
            "homePage": "https://github.com/andrecossa5/nf-lenti.git"
        },
        "run_info": {
            "date": datetime.datetime.now().strftime("%Y-%m-%d"),
            "user": user_name,
            "working_directory": working_dir,
        },
        "parameters": {
            "indir": args.indir,
            "outdir": args.params_outdir,
            "anchor_sequence": args.anchor_sequence,
            "min_n_reads": args.min_n_reads,
            "hamming_treshold": args.hamming_treshold,
        },
        "sample": args.sample,
        "metrics": {
            "total_GBC_reads": int(total_GBC_reads),
            "n_unique_GBCs": int(n_unique_GBCs),
            "n_filtered_GBCs": int(n_filtered_GBCs),
            "n_corrected_GBCs": int(n_corrected),
            "n_degenerated_GBCs": int(n_degenerated),
            "median_n_reads_added": float(median_n_reads_added),
        },
        "metrics_info": {
            "total_GBC_reads": "Total reads assigned to valid GBCs (18 bp) before filtering and correction.",
            "n_unique_GBCs": "Unique GBC (18bp) detected before filtering and correction.",
            "n_filtered_GBCs": "Final GBCs retained after spike-in removal, correction, and read threshold filtering.",
            "n_corrected_GBCs": "GBCs kept as correction targets after collapsing sequences based on Hamming distance clustering. Default: 3",
            "n_degenerated_GBCs": "GBCs merged into corrected barcodes.",
            "median_n_reads_added": "Median reads reassigned to corrected GBCs from degenerate barcodes."
        },
        "qc_images": qc_encoded
    }

    # Save JSON
    json_path_outdir = Path(args.outdir) / 'run_summary.json'
    with open(json_path_outdir, 'w') as jf:
        json.dump(summary_data, jf, indent=2)

if __name__ == '__main__':
    main()
