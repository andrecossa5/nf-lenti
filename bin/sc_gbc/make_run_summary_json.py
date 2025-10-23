#!/usr/bin/env python3
"""
make_run_summary_json.py

Generate run_summary.json for nf-lenti single-cell GBC pipeline.
- Collects metrics from clone calling summary, cells summary, clones, GBCs, and bulk reference files.
- Optionally counts reads from FASTQ or BAM input.
- Encodes QC images for embedding in the report.
- Outputs a JSON summary with pipeline, run info, parameters, metrics, and QC images.
"""

import argparse, json, base64, os, pandas as pd
from datetime import date
from pathlib import Path

import sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "common"))
from utils import (
    encode_image_b64,
    resolve_working_dir_and_user,
    line_count,
    count_unique_gbc_in_sc,
    count_fastq_reads,
    count_bam_reads,
    parse_clone_calling_summary,
)

def read_html_content(html_path):
    """Read HTML file content as a string, or return None if not found."""
    if html_path and Path(html_path).exists():
        with open(html_path, "r", encoding="utf-8") as f:
            return f.read()
    return None

def encode_html_b64(html_path):
    """Base64-encode HTML file content, or return None if not found."""
    if html_path and Path(html_path).exists():
        with open(html_path, "rb") as f:
            return base64.b64encode(f.read()).decode("utf-8")
    return None

def get_bam_path_from_csv(bam_csv, sample):
    df = pd.read_csv(bam_csv)
    row = df.loc[df["sample"] == sample]
    return row["bam"].iloc[0] if not row.empty else None

def main():
    # Create the parser
    parser = argparse.ArgumentParser(
        prog="make_run_summary_json",
        description="""
        Generate run_summary.json for nf-lenti single-cell GBC pipeline.
        - Collects metrics from clone calling summary, cells summary, clones, GBCs, and bulk reference files.
        - Optionally counts reads from FASTQ or BAM input.
        - Encodes QC images for embedding in the report.
        - Outputs a JSON summary with pipeline, run info, parameters, metrics, and QC images.
        """
    )
    parser.add_argument(
        "--sample",
        type=str,
        required=True,
        help="Sample name."
    )
    parser.add_argument(
        "--clone_summary_txt",
        type=str,
        required=True,
        help="Path to clone_calling_summary.txt."
    )
    parser.add_argument(
        "--cells_summary",
        type=str,
        required=True,
        help="Path to cells_summary file."
    )
    parser.add_argument(
        "--clones",
        type=str,
        required=True,
        help="Path to clones file."
    )
    parser.add_argument(
        "--gbcs",
        type=str,
        required=True,
        help="Path to GBCs file."
    )
    parser.add_argument(
        "--bulk_gbc",
        type=str,
        required=True,
        help="Path to bulk GBC reference file."
    )
    parser.add_argument(
        "--combo_plot",
        type=str,
        help="Path to combo plot image."
    )
    parser.add_argument(
        "--umi_dist",
        type=str,
        help="Path to UMI distribution image."
    )
    parser.add_argument(
        "--moi_dist",
        type=str,
        help="Path to MOI distribution image."
    )
    parser.add_argument(
        "--clone_sz",
        type=str,
        help="Path to clone size image."
    )
    parser.add_argument(
        "--umi_dist_interactive",
        type=str,
        help="Path to interactive UMI distribution image."
    )
    parser.add_argument(
        "--moi_dist_interactive",
        type=str,
        help="Path to interactive MOI distribution image."
    )
    parser.add_argument(
        "--clone_sz_interactive",
        type=str,
        help="Path to interactive clone size image."
    )
    parser.add_argument(
        "--raw_data_input",
        type=str,
        help="Path to raw data input folder or file."
    )
    parser.add_argument(
        "--raw_data_input_type",
        type=str,
        help="Type of raw data input (fastq, fastq,gbc, bam)."
    )
    parser.add_argument(
        "--sc_outdir",
        type=str,
        help="Path to single-cell output directory."
    )
    parser.add_argument(
        "--pattern",
        type=str,
        help="Pattern used for demultiplexing."
    )
    parser.add_argument(
        "--ref",
        type=str,
        help="Reference used."
    )
    parser.add_argument(
        "--out_json",
        type=str,
        default="run_summary.json",
        help="Output JSON file name. Default: run_summary.json."
    )

    args = parser.parse_args()

    # --- Metrics ---
    metrics = parse_clone_calling_summary(args.clone_summary_txt)
    metrics.setdefault("unique_good_GBCs", None)
    metrics["unique_gbc_bulk_ref"] = line_count(args.bulk_gbc)
    metrics["unique_gbc_sc"] = count_unique_gbc_in_sc(args.gbcs)

    # Add total reads
    if args.raw_data_input and args.raw_data_input_type:
        input_type = args.raw_data_input_type.lower()
        if input_type in ["fastq", "fastq,gbc"]:
            metrics["total_reads"] = count_fastq_reads(args.raw_data_input)
        elif input_type == "bam":
            bam_path = args.raw_data_input
            if bam_path.endswith(".csv"):
                bam_path = get_bam_path_from_csv(bam_path, args.sample)
            if bam_path:
                mapped, unmapped = count_bam_reads(bam_path)
                metrics["mapped_reads"] = mapped
                metrics["unmapped_reads"] = unmapped

    # STARsolo outputs
    if args.sc_outdir:
        barcodes = Path(args.sc_outdir) / "filtered" / "barcodes.tsv.gz"
        if barcodes.exists():
            metrics["n_putative_cells"] = line_count(barcodes)
        transcripts = Path(args.sc_outdir) / "filtered" / "features.tsv.gz"
        if transcripts.exists():
            metrics["total_n_transcripts"] = line_count(transcripts)

    working_dir, user_name = resolve_working_dir_and_user()

    # QC images as base64
    qc_image_info = [
        (args.combo_plot, "CBC-GBC Combo Support"),
        (args.umi_dist, "UMI Count Distribution per CBC-GBC"),
        (args.moi_dist, "MOI Distribution"),
        (args.clone_sz, "Clone Size Distribution"),
    ]
    interactive_map = {
        "UMI Count Distribution per CBC-GBC": args.umi_dist_interactive,
        "MOI Distribution": args.moi_dist_interactive,
        "Clone Size Distribution": args.clone_sz_interactive
    }
    qc_images = []
    for pth, title in qc_image_info:
        if pth:
            enc = encode_image_b64(pth)
            if enc:
                img_dict = {"image": enc, "title": title}
                interactive_html_path = interactive_map.get(title)
                if interactive_html_path:
                    html_b64 = encode_html_b64(interactive_html_path)
                    if html_b64:
                        img_dict["interactive_html_b64"] = html_b64
                qc_images.append(img_dict)

    report = {
        "run_type": "sc_gbc",
        "pipeline": {
            "name": "nf-lenti",
            "version": "1.0.0",
            "homePage": "https://github.com/andrecossa5/nf-lenti.git"
        },
        "run_info": {
            "date": date.today().isoformat(),
            "user": user_name,
            "working_directory": working_dir,
        },
        "parameters": {
            "raw_data_input": args.raw_data_input,
            "raw_data_input_type": args.raw_data_input_type,
            "sc_outdir": args.sc_outdir,
            "pattern": args.pattern,
            "ref": args.ref,
        },
        "metrics_info": {
            "total_reads": "Total number of sequencing reads",
            "mapped_reads": "Number of mapped reads in BAM",
            "unmapped_reads": "Number of unmapped reads in BAM",
            "unique_gbc_bulk_ref": "Unique GBC in bulk reference",
            "unique_gbc_sc": "Unique GBC detected in single-cell data",
            "unique_good_GBCs": "Unique GBCs retained after quality filtering",
            "unsupported_combos": "CBC-GBC combos unsupported",
            "supported_combos": "CBC-GBC combos supported",
            "observed_moi_median": "Median observed MOI",
            "observed_moi_std": "Std of observed MOI",
            "starting_CBCs": "Number of starting CBCs",
            "uniquely_barcoded_cells": "Number of cells that are confidently assigned to a single GBC",
            "n_clones": "The number of unique GBCs found among cells that have only one GBC",
            "n_clones_10plus": "Number of clones (unique GBCs) that have 10 or more cells",
            "n_putative_cells": "Putative cells from STARsolo",
            "total_n_transcripts": "Number of transcripts detected (STARsolo)"
        },
        "samples": [
            {
                "id": args.sample,
                "metrics": metrics,
                "qc_images": qc_images,
                "raw_data_input_type": args.raw_data_input_type
            }
        ]
    }
    with open(args.out_json, "w") as f:
        json.dump(report, f, indent=2)

if __name__ == "__main__":
    main()