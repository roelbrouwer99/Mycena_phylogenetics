# =========================
# README
# =========================
# Purpose: Create a tanglegram between two unrooted trees (SVG output).
# Inputs:
#   - tree1 and tree2 file paths (tree1, tree2)
# Outputs:
#   - SVG file (path in svg())
# How to run:
#   Rscript Tanglegram_unrooted.R
# Key settings to tweak:
#   tree1, tree2, file1, file2, output SVG path.
# Notes:
#   - diff_labels should be defined to color differing tips.
# =========================

library(ape)
library(phytools)

# --- Load rooted trees ---
tree1 < -read.tree("C:/Users/rbrou/Desktop/Roel/Naturalis/Phylogenetics/3_Tree files/BOLD/Tanglegrams/BOLD_mycenaceae_FFT-NS-2_untrimmed.treefile")
tree2 < -read.tree("C:/Users/rbrou/Desktop/Roel/Naturalis/Phylogenetics/3_Tree files/BOLD/Tanglegrams/BOLD_mycenaceae_FFT-NS-2_manual_ends.treefile")

# --- Input file names ---
file1 < -"BOLD_mycenaceae_FFT-NS-2_untrimmed.treefile"
file2 < -"BOLD_mycenaceae_FFT-NS-2_manual_ends.treefile"

# --- Rooted trees already aligned by outgroup
# Ensure identical tip order
tree2 < -keep.tip(tree2, tree1$tip.label)
tree2 < -reorder.phylo(tree2, order(match(tree2$tip.label, tree1$tip.label)))

# --- Use cophylo to align trees ---
cophy < -cophylo(tree1, tree2, rotate = TRUE)

# --- Ensure linking lines vector is correct length ---
link_colors < -rep("grey40", Ntip(tree1)) # default color

# Tip label colors
tip_colors_left < -ifelse(cophy$tree1$tip.label % in % diff_labels, "red", "black")
tip_colors_right < -ifelse(cophy$tree2$tip.label % in % diff_labels, "red", "black")

# --- Save SVG instead of PDF ---
svg("C:/Users/rbrou/Desktop/Roel/Naturalis/Phylogenetics/3_Tree files/BOLD/Tanglegrams/Untrimmed_Manual.svg", width = 35, height = 30)

plot(cophy,
     link.type = "curved",
     link.lwd = 0.5,
     link.col = link_colors,
     fsize = 0.3
)

# Add colored tip labels
tiplabels(pch = NA, col = tip_colors_left, tree = cophy$tree1)
tiplabels(pch = NA, col = tip_colors_right, tree = cophy$tree2)

# Add titles
mtext(side = 3, line = 2, at = 0.25, text = paste0("Tree 1: ", file1), cex = 1.2)
mtext(side = 3, line = 2, at = 0.75, text = paste0("Tree 2: ", file2), cex = 1.2)

dev.off()
cat("✅ Tanglegram with connecting lines saved as SVG\n")
