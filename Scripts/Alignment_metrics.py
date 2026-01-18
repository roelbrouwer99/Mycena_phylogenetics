#!/usr/bin/env python3
# =========================
# README
# =========================
# Purpose: Composite figure summarizing alignment and tree quality statistics.
# Inputs:
#   - alignment_tree_metrics.csv (alignment/tree metrics)
#   - bootstrap_summary.csv (bootstrap summaries)
# Outputs:
#   - composite_alignment_tree_figure.png
#   - PanelE_bootstrap_violin.png
#   - PanelF_heatmap_summary.png
#   - All outputs saved to the output_dir folder.
# How to run:
#   python alignment_metrics.py
# Key settings to tweak:
#   output_dir, input CSV paths, figure sizes, and scale_factor.
# Notes:
#   - Requires pandas, numpy, matplotlib, seaborn, and biopython.
# =========================


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from Bio import Phylo

# ----------------------------------------------------------
# 1. Create output folder
# ----------------------------------------------------------
output_dir = Path("C:/Users/rbrou/Desktop/Alignment_metrics")
output_dir.mkdir(exist_ok=True)
print(f"Output folder created or already exists: {output_dir.resolve()}")

# ----------------------------------------------------------
# 2. Load data
# ----------------------------------------------------------
# File 1: alignment/tree metrics
metrics = pd.read_csv(
    "C:/Users/rbrou/Desktop/Naturalis scripts/alignment_tree_metrics.csv", sep=";", decimal=",")
# Expected columns: method, length, RCV, evol_rate, treeness, treeness_RCV, mean_boot, [trimming]

metrics.columns = (
    metrics.columns.str.strip()
    .str.lower()
    .str.replace(" ", "_")
)

metrics.rename(columns={
    "alignment_method": "method",
    "var._sites": "var_sites",
    "internal_branch_statistics_(mean)": "ibs_mean",
    "terminal_branch_statistics_(mean)": "tbs_mean",
    "ibs/tbs": "ibs_tbs",
    "treeness_over_rcv": "treeness_rcv"
}, inplace=True)

print("✅ Columns loaded:", list(metrics.columns))

# File 2: bootstrap summaries
boots = pd.read_csv(
    "C:/Users/rbrou/Desktop/Naturalis scripts/bootstrap_summary.csv", sep=";", decimal=",")

# ----------------------------------------------------------
# 3. Panel A–C: Alignment quality scatterplots
# ----------------------------------------------------------
sns.set(style="whitegrid", context="talk")

fig, axes = plt.subplots(2, 2, figsize=(13, 10))

# Panel A – Overall trend
axA = axes[0, 0]
sns.scatterplot(data=metrics, x="length", y="rcv", ax=axA, color="black", s=40)
sns.regplot(data=metrics, x="length", y="rcv", scatter=False, ax=axA,
            color="black", line_kws={'linestyle': 'dotted'})
axA.set_title("A  Alignment length vs. RCV (overall)")
axA.set_xlabel("Alignment length")
axA.set_ylabel("RCV score")

# Panel B – Color by alignment method
axB = axes[0, 1]
sns.scatterplot(data=metrics, x="length", y="rcv", hue="method", ax=axB, s=60)
axB.legend(bbox_to_anchor=(1.05, 1), title="Alignment method")
axB.set_title("B  Colored by alignment method")
axB.set_xlabel("Alignment length")
axB.set_ylabel("RCV score")

# Panel C – Color by trimming (if present)
axC = axes[1, 0]
if "trimming" in metrics.columns:
    sns.scatterplot(data=metrics, x="length", y="rcv", hue="trimming",
                    style="trimming", ax=axC, s=70)
    axC.legend(bbox_to_anchor=(1.05, 1), title="Trimming")
    axC.set_title("C  Colored by trimming strategy")
else:
    axC.text(0.5, 0.5, "No trimming data provided", ha="center", va="center",
             fontsize=12, color="gray", transform=axC.transAxes)
axC.set_xlabel("Alignment length")
axC.set_ylabel("RCV score")

# ----------------------------------------------------------
# 4. Panel D: Radar plot of tree signal metrics
# ----------------------------------------------------------
df_norm = metrics.copy()
cols = ["treeness", "treeness_rcv", "evolutionary_rate", "rcv"]

# normalize columns
df_norm[cols] = df_norm[cols].apply(
    lambda x: (x - x.min()) / (x.max() - x.min()))

# invert measures where lower = better
df_norm["rcv"] = 1 - df_norm["rcv"]
df_norm["evolutionary_rate"] = 1 - df_norm["evolutionary_rate"]


labels = cols
angles = np.linspace(0, 2*np.pi, len(labels), endpoint=False).tolist()
angles += angles[:1]

axR = axes[1, 1]
axR = plt.subplot(2, 2, 4, polar=True)
for _, row in df_norm.iterrows():
    vals = row[labels].tolist()
    vals += vals[:1]
    axR.plot(angles, vals, label=row["method"], linewidth=1)
    axR.fill(angles, vals, alpha=0.1)
axR.set_xticks(angles[:-1])
axR.set_xticklabels(labels)
axR.set_yticklabels([])
axR.set_title("D  Normalized tree-signal metrics (Radar)")
plt.legend(bbox_to_anchor=(1.4, 1.1))

plt.tight_layout()
plt.savefig(output_dir / "composite_alignment_tree_figure.png", dpi=300)
print(
    f"Saved composite figure: {output_dir / 'composite_alignment_tree_figure.png'}")

# ----------------------------------------------------------
# 5. Panel E: Violin plot of bootstrap supports
# ----------------------------------------------------------
# Prepare data
boots_long = pd.read_csv("bootstrap_summary.csv")
boots_long.columns = [c.strip().lower() for c in boots_long.columns]

if "method" in boots_long.columns:
    fig, ax = plt.subplots(figsize=(10, 5))
    sns.violinplot(data=boots_long, x="method", y="mean",
                   inner="box", palette="Set2", cut=0)
    plt.ylim(0, 100)
    plt.axhline(80, ls="--", color="gray", lw=1)
    plt.axhline(95, ls="--", color="gray", lw=1)
    plt.ylabel("Bootstrap support (%)")
    plt.title("E  Bootstrap support distributions by alignment method (Violin)")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.savefig(output_dir / "PanelE_bootstrap_violin.png", dpi=300)
    print(
        f"Saved Panel E violin plot: {output_dir / 'PanelE_bootstrap_violin.png'}")
else:
    print("⚠️ Skipped Panel E: bootstrap_summary.csv missing 'method' column")

# ----------------------------------------------------------
# 6. Heatmap Summary
# ----------------------------------------------------------
heat = df_norm.set_index("method")[cols]
plt.figure(figsize=(6, 4))
sns.heatmap(heat, cmap="RdYlGn", annot=True, cbar=False)
plt.title("Summary ranking (higher = better)")
plt.tight_layout()
plt.savefig(output_dir / "PanelF_heatmap_summary.png", dpi=300)
print(f"Saved heatmap: {output_dir / 'PanelF_heatmap_summary.png'}")

print("\nAll outputs saved to:", output_dir.resolve())
