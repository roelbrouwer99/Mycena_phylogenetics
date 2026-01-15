#!/usr/bin/env python3
# =========================
# README
# =========================
# Purpose: Summarize bootstrap support values from IQ-TREE .contree files.
# Inputs:
#   - One or more IQ-TREE .contree files passed on the command line.
# Outputs:
#   - bootstrap_summary.csv (mean, median, prop80, prop95, n_nodes)
# How to run:
#   python Bootstrap_analysis.py *.contree
# Key settings to tweak:
#   None (the script summarizes whatever files you pass in).
# Notes:
#   - Requires numpy, pandas, and biopython.
# =========================


import sys
import numpy as np
import pandas as pd
from Bio import Phylo
from pathlib import Path


def parse_bootstrap_supports(tree_path):
    """Extract bootstrap values from an IQ-TREE .contree file."""
    tree = Phylo.read(tree_path, "newick")
    supports = []
    for clade in tree.find_clades():
        if clade.confidence is not None:
            supports.append(clade.confidence)
    return np.array(supports)


def summarize_supports(values):
    """Return summary statistics for a bootstrap vector."""
    if len(values) == 0:
        return dict(mean=np.nan, median=np.nan, prop80=np.nan, prop95=np.nan)
    return dict(
        mean=np.mean(values),
        median=np.median(values),
        prop80=np.mean(values >= 80),
        prop95=np.mean(values >= 95)
    )


def main(paths):
    summary_rows = []
    for path in paths:
        vals = parse_bootstrap_supports(path)
        stats = summarize_supports(vals)
        stats["method"] = Path(path).stem
        stats["n_nodes"] = len(vals)
        summary_rows.append(stats)
        print(f"{Path(path).name}: mean={stats['mean']:.2f}, "
              f"median={stats['median']:.1f}, ≥80%={stats['prop80']:.2%}, "
              f"≥95%={stats['prop95']:.2%}")

    df = pd.DataFrame(summary_rows).set_index(
        "method").sort_values("mean", ascending=False)
    df.to_csv("bootstrap_summary.csv", float_format="%.3f")
    print("\nSaved summary table: bootstrap_summary.csv")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit("Usage: python analyze_bootstrap_supports_tableonly.py *.contree")
    paths = sys.argv[1:]
    main(paths)
