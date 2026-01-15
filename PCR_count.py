# =========================
# README
# =========================
# Purpose: Plot PCR success/failure counts over collection years with stockplate bands.
# Inputs:
#   - PCR_count.csv (expects PCR_condition, year, count)
# Outputs:
#   - Interactive plot shown on screen (no file saved by default).
# How to run:
#   python PCR_count.py
# Key settings to tweak:
#   CSV path, stockplate_groups, pcr_true_names, x-axis limits, scale_factor.
# Notes:
#   - Requires pandas, numpy, and matplotlib.
# =========================

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch

plt.style.use('seaborn-v0_8-whitegrid')

# Load CSV
df = pd.read_csv("C:/Users/rbrou/Desktop/PCR_count/PCR_count.csv", sep=";")

# Extract dataset base names (removing Worked/Failed) and reverse sort
datasets = sorted(set([cond.replace(' - Failed', '').replace(' - Worked!', '')
                       for cond in df["PCR_condition"]]), reverse=True)
dataset_pcr = {dataset: ' '.join(
    dataset.split(' ')[:2]) for dataset in datasets}

# Map PCR numbers to true names
pcr_true_names = {
    "PCR 1": "LSU (long)",
    "PCR 2": "LSU (short)",
    "PCR 3": "ITS",
    "PCR 4": "LSU (long)",
    "PCR 5": "LSU (short)",
    "PCR 6": "ITS",
    "PCR 7": "LSU (long)",
    "PCR 8": "LSU (short)",
    "PCR 9": "ITS",
    "PCR 10": "LSU (long)",
    "PCR 11": "LSU (short)",
    "PCR 12": "ITS"
}

# Define which PCRs belong to which stockplate
stockplate_groups = {
    "NCBN002483": ["PCR 1", "PCR 2", "PCR 3"],
    "NCBN002846": ["PCR 4", "PCR 5", "PCR 6"],
    "NCBN002950": ["PCR 7", "PCR 8", "PCR 9"],
    "NCBN002951": ["PCR 10", "PCR 11", "PCR 12"]
}

colors = {'Worked': 'forestgreen', 'Failed': 'firebrick'}
legend_handles = {}

fig, ax = plt.subplots(figsize=(12, 7))

scale_factor = 0.35
gap = 1.5  # fixed vertical gap between PCRs

# --- Compute cumulative y-positions to prevent spike overlaps ---
y_positions = {}
y_top = 0  # starting baseline at the top

for dataset in datasets:
    sub_worked = df[df["PCR_condition"] == dataset + ' - Worked!']
    sub_failed = df[df["PCR_condition"] == dataset + ' - Failed']

    max_up = sub_worked['count'].max(
    ) * scale_factor if not sub_worked.empty else 0
    max_down = sub_failed['count'].max(
    ) * scale_factor if not sub_failed.empty else 0

    # Baseline for this PCR
    y_positions[dataset] = y_top - max_up - gap

    # Update y_top for next PCR to include downward spike
    y_top = y_positions[dataset] - max_down

# --- Draw stockplate bands with rounded rectangles ---
padding = 0.8  # vertical padding above/below PCRs
rounding = 3   # corner rounding size

# Manual x-coordinates
x_left = 1875
x_right = 2045

for sp_name, pcrs in stockplate_groups.items():
    related = [d for d in datasets if any(d.startswith(p) for p in pcrs)]
    if not related:
        continue

    y_top_band = y_positions[related[0]] + padding
    y_bottom_band = y_positions[related[-1]] - padding
    height = y_top_band - y_bottom_band

    rect = FancyBboxPatch(
        (x_left, y_bottom_band),
        x_right - x_left, height,
        boxstyle=f"round,pad=0.02,rounding_size={rounding}",
        linewidth=0, facecolor='lightgray', alpha=0.25, zorder=0
    )
    ax.add_patch(rect)

    # Vertical stockplate label
    ax.text(x_left + 5, (y_top_band + y_bottom_band) / 2, sp_name,
            ha='left', va='center', fontsize=10, fontweight='bold',
            color='black', rotation=90)

# --- Plot PCR datasets ---
for i, dataset in enumerate(datasets):
    ypos = y_positions[dataset]
    pcr_group = dataset_pcr[dataset]

    sub_worked = df[df["PCR_condition"] == dataset + ' - Worked!']
    sub_failed = df[df["PCR_condition"] == dataset + ' - Failed']

    all_years = pd.concat([sub_worked['year'], sub_failed['year']])
    start_year = all_years.min()
    end_year = all_years.max()

    # Alternating gray timeline bars
    bar_color = "#c5c5c5" if i % 2 == 0 else "#959595"
    ax.plot([start_year, end_year], [ypos, ypos], lw=5,
            solid_capstyle='round', color=bar_color, alpha=1)

    # Spikes above (worked)
    for _, row in sub_worked.iterrows():
        ax.vlines(x=row['year'], ymin=ypos, ymax=ypos + row['count'] * scale_factor,
                  color=colors['Worked'], lw=3, alpha=0.9)

    # Spikes below (failed)
    for _, row in sub_failed.iterrows():
        ax.vlines(x=row['year'], ymin=ypos - row['count'] * scale_factor, ymax=ypos,
                  color=colors['Failed'], lw=3, alpha=0.9)

    # Success rate (%)
    total_passed = sub_worked['count'].sum()
    total_failed = sub_failed['count'].sum()
    total = total_passed + total_failed
    success_rate = (total_passed / total * 100) if total > 0 else np.nan

    # Label: PCR true name + percentage below
    true_name = pcr_true_names.get(pcr_group, pcr_group)
    label_text = f"{true_name}\n({success_rate:.1f}%)"  # divide over two lines
    ax.text(end_year + 2, ypos, label_text,
            ha='left', va='center', fontsize=10, fontweight='bold', color='black')

    # Legend handles
    if 'Worked' not in legend_handles and not sub_worked.empty:
        legend_handles['Worked'] = ax.plot([], [], color=colors['Worked'], lw=3,
                                           label='Passed samples')[0]
    if 'Failed' not in legend_handles and not sub_failed.empty:
        legend_handles['Failed'] = ax.plot([], [], color=colors['Failed'], lw=3,
                                           label='Failed samples')[0]

# --- Scale bar ---
ref_y = y_top - 3
ref_x = 1875
ref_height = 5 * scale_factor
ax.vlines(x=ref_x, ymin=ref_y, ymax=ref_y + ref_height, color='black', lw=2)
ax.text(ref_x + 2, ref_y + ref_height / 2,
        "= 5 samples", va='center', fontsize=9, color='black', fontweight='bold')

# --- Formatting ---
ax.set_yticks([])
ax.set_xlabel("Collection Year", fontsize=11, color='black', fontweight='bold')
ax.set_xlim(1870, 2050)
ax.set_xticks(range(1870, 2051, 10))
ax.set_title("Passed vs Failed Samples based on collection year",
             fontsize=13, pad=15, color='black', fontweight='bold')
ax.grid(axis='x', linestyle='--', linewidth=0.5, alpha=0.6)

# Legend (upper left)
ax.legend(handles=list(legend_handles.values()), loc='upper left',
          frameon=True, fontsize=9)

plt.tight_layout()
plt.show()
