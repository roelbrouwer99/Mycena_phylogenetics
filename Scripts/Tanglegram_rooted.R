# =========================
# README
# =========================
# Purpose: Create a rooted ML vs BI tanglegram with outgroup rooting.
# Inputs:
#   - ML tree file (tree_path_A)
#   - BI tree file (tree_path_B)
#   - Outgroup taxa labels (outgroup_taxa)
# Outputs:
#   - PDF tanglegram (output_pdf)
# How to run:
#   Rscript Tanglegram_rooted.R
# Key settings to tweak:
#   tree_path_A, tree_path_B, output_pdf, outgroup_taxa.
# Notes:
#   - tip_col_A and tip_col_B must be defined if colored tips are desired.
# =========================

library(ape)
library(phytools)
library(phangorn)

#------------------------------------------------------------
# File paths
#------------------------------------------------------------
tree_path_A <- "C:/Users/rbrou/Desktop/Roel/Naturalis/Phylogenetics/3_Tree files/PRANK_Trees/IQ-TREE/contree/Perry_ITS_PRANK.best.fas.clipkit.contree"
tree_path_B <- "C:/Users/rbrou/Documents/cluster_3_consensus.tre"
output_pdf  <- "C:/Users/rbrou/Desktop/Roel/Naturalis/Phylogenetics/4_Tanglegrams/BI_clusters/ML_BIcl3_rooted.pdf"

#------------------------------------------------------------
# Load trees
#------------------------------------------------------------
TreeA <- read.tree(tree_path_A)  # ML
TreeB <- read.tree(tree_path_B)  # BI

#------------------------------------------------------------
# Root using exact outgroup labels
#------------------------------------------------------------
outgroup_taxa <- c(
  "NR_184518_Effusomyces_thailandicus_BJFC_023496_Thailand_Type",
  "NR_184519_Parvodontia_austrosinensis_BJFC_024249_China_Type",
  "ON117192_Cystostereum_submurrayi_He_4379_China",
  "NR_184516_Crustomyces_albidus_BJFC_033109_China_Type",
  "DQ486687_Nothopanus_candidissimus_PBM_2411__WTU__Washington",
  "MH857243_Chondrostereum_purpureum_CBS_352_53_France",
  "KY614001_Gloeostereum_incarnatum_BCC_41461_Thailand",
  "AY230866_Campanophyllum_proboscideum_TENN56402_Mexico",
  "DQ486698_Cyphella_digitalis_CBS679_82"
)

outgroupA <- intersect(outgroup_taxa, TreeA$tip.label)
outgroupB <- intersect(outgroup_taxa, TreeB$tip.label)

TreeA_rooted <- root(TreeA, outgroup = outgroupA, resolve.root = TRUE)
TreeB_rooted <- root(TreeB, outgroup = outgroupB, resolve.root = TRUE)

#------------------------------------------------------------
# Orientation: root at bottom
#------------------------------------------------------------
TreeA_rooted <- ladderize(TreeA_rooted, right = FALSE)
TreeB_rooted <- ladderize(TreeB_rooted, right = FALSE)

#------------------------------------------------------------
# Clean labels for association matching
#------------------------------------------------------------
clean_label <- function(x) {
  x <- gsub("[^A-Za-z0-9]", "", x)   # keep alphanumerics
  toupper(x)
}

tipsA_clean <- clean_label(TreeA_rooted$tip.label)
tipsB_clean <- clean_label(TreeB_rooted$tip.label)

shared_clean <- intersect(tipsA_clean, tipsB_clean)

if (length(shared_clean) == 0)
  stop("No shared taxa found after label cleaning.")

idxA <- match(shared_clean, tipsA_clean)
idxB <- match(shared_clean, tipsB_clean)

association <- cbind(
  TreeA_rooted$tip.label[idxA],
  TreeB_rooted$tip.label[idxB]
)

#------------------------------------------------------------
# Build cophylogeny object
#------------------------------------------------------------
obj <- cophylo(TreeA_rooted, TreeB_rooted,
               assoc = association,
               rotate = FALSE,
               plot = FALSE)

#------------------------------------------------------------
# Link colors (movement)
#------------------------------------------------------------
left_order  <- obj$trees[[1]]$tip.label
right_order <- obj$trees[[2]]$tip.label
rank_left   <- match(obj$assoc[,1], left_order)
rank_right  <- match(obj$assoc[,2], right_order)
rank_diff   <- abs(rank_left - rank_right)

palette_fun <- colorRampPalette(c("forestgreen", "gold", "firebrick"))
n_col <- 200

link_colors <- palette_fun(n_col)[
  round(rank_diff / max(rank_diff) * (n_col - 1)) + 1
]

#------------------------------------------------------------
# PLOT EVERYTHING
#------------------------------------------------------------
pdf(output_pdf, height = 24, width = 18)
par(mar = c(5, 5, 5, 2))

plot(obj,
     link.type       = "curved",
     link.lwd        = 3,
     link.lty        = "solid",
     link.col        = make.transparent(link_colors, 0.8),
     fsize           = 0.5,
     label.offset    = 0,
     tip.color.left  = tip_col_A,
     tip.color.right = tip_col_B)


title("ML vs BI-cl3",
      font.main = 3)

dev.off()

cat("Plot saved to:/n", output_pdf, "/n")


