"""
utils.py

Utility functions for nf-lenti pipelines.
- Image encoding for HTML/JSON embedding.
- Working directory and user resolution.
- Metrics extraction from summary files.
- File line counting and barcode counting utilities.
- FASTQ/BAM read counting.
- Parsing clone calling summaries.
"""

import base64
from pathlib import Path
import os
from io import BytesIO
import matplotlib.pyplot as plt

def encode_image_b64(path_img):
    p = Path(path_img)
    with open(p, "rb") as f:
        b64 = base64.b64encode(f.read()).decode("utf-8")
    ext = p.suffix.lower().lstrip(".")
    if ext == "jpg":
        ext = "jpeg"
    return f"data:image/{ext};base64,{b64}"

def resolve_working_dir_and_user():
    task_dir = Path(os.getcwd()).resolve()
    working_dir = next((str(parent.parent) for parent in task_dir.parents if parent.name == "work"), str(task_dir))
    user_name = os.environ.get("USER") or "unknown"
    parts = Path(working_dir).parts
    if "Users" in parts:
        idx = parts.index("Users")
        if idx + 1 < len(parts):
            candidate = parts[idx + 1]
            if candidate not in ("work", "nf-core", "root") and not candidate.startswith("conda"):
                user_name = candidate
    return working_dir, user_name

def encode_plot_to_base64(fig):
    """Convert Matplotlib figure to Base64-encoded PNG for embedding in JSON/HTML."""
    buf = BytesIO()
    fig.savefig(buf, format="png", dpi=150, bbox_inches="tight")
    buf.seek(0)
    b64 = base64.b64encode(buf.read()).decode("utf-8")
    plt.close(fig)
    return f"data:image/png;base64,{b64}"

def get_sample_id_from_json(data, fallback):
    """Extract a single sample ID from run_summary.json, else fallback to folder name."""
    samples = data.get("samples") or []
    if len(samples) == 1 and samples[0].get("id"):
        return samples[0]["id"]
    for s in samples:
        if s.get("id"):
            return s["id"]
    return fallback

def get_metrics_for_sample(data, sample_id):
    """Return metrics dict for the matching sample_id"""
    samples = data.get("samples") or []
    for s in samples:
        if s.get("id") == sample_id:
            return s.get("metrics") or {}
    for s in samples:
        m = s.get("metrics")
        if m:
            return m
    return {}

def line_count(pth):
    if not pth: return None
    p = Path(pth)
    if not p.exists(): return None
    n = 0
    if p.suffix == ".gz":
        import gzip
        with gzip.open(p, "rt", newline="") as f:
            for _ in f: n += 1
    else:
        with open(p, "r", newline="") as f:
            for _ in f: n += 1
    return n

def count_unique_gbc_in_sc(path):
    import pandas as pd
    df = pd.read_csv(path, sep="\t", compression="gzip")
    return df["GBC"].nunique() if "GBC" in df.columns else None

def count_fastq_reads(fastq_folder):
    import gzip
    total_reads = 0
    for fname in os.listdir(fastq_folder):
        if "R1" in fname and fname.endswith(".fastq.gz"):
            with gzip.open(os.path.join(fastq_folder, fname), "rt") as fh:
                total_reads += sum(1 for _ in fh) // 4
    return total_reads

def count_bam_reads(bam_path):
    import pysam
    bam = pysam.AlignmentFile(bam_path, "rb")
    mapped = unmapped = 0
    for read in bam:
        if read.is_unmapped:
            unmapped += 1
        else:
            mapped += 1
    bam.close()
    return mapped, unmapped

def parse_clone_calling_summary(path):
    import re
    metrics = {}
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith("- Unique, \"good\" GBCs:"):
                metrics["unique_good_GBCs"] = int(re.findall(r"\d+", line)[0])
            elif line.startswith("- n unsupported CBC-GBC combos"):
                metrics["unsupported_combos"] = int(re.findall(r"\d+", line)[0])
            elif line.startswith("- n supported CBC-GBC combos"):
                metrics["supported_combos"] = int(re.findall(r"\d+", line)[0])
            elif line.startswith("- Observed MOI"):
                match = re.search(r"([\d\.]+) median, \(?\+\-([\d\.]+)\)?", line)
                if match:
                    metrics["observed_moi_median"] = float(match.group(1))
                    metrics["observed_moi_std"] = float(match.group(2))
            elif line.startswith("- n starting CBC"):
                metrics["starting_CBCs"] = int(re.findall(r"\d+", line)[0])
            elif line.startswith("- n uniquely barcoded cells"):
                metrics["uniquely_barcoded_cells"] = int(re.findall(r"\d+", line)[0])
            elif line.startswith("- n clones:"):
                metrics["n_clones"] = int(re.findall(r"\d+", line)[0])
            elif line.startswith("- n clones>=10 cells"):
                matches = re.findall(r"\d+", line)
                if matches:
                    metrics["n_clones_10plus"] = int(matches[-1])
    return metrics
