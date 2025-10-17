#!/usr/bin/env Rscript

# -----------------------------
# Load required libraries
# -----------------------------
library(TreeDist)
library(ape)

# -----------------------------
# Command-line arguments
# -----------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript compute_spr_distance.R <input_newick_file> <output_csv>")
}

input_file <- args[1]
output_file <- args[2]

# -----------------------------
# Read input file
# -----------------------------
# Expected format per line:
# start <tab> end <tab> newick_string
trees_df <- read.table(input_file, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(trees_df) <- c("start", "end", "newick")

# -----------------------------
# Compute SPR distance between adjacent trees
# -----------------------------
results <- data.frame(start = numeric(), end = numeric(), spr_dist = numeric(),
                      splitInfoDist = numeric(), rf_dist = numeric(), stringsAsFactors = FALSE

for (i in 1:(nrow(trees_df) - 1)) {
  line1_newick <- trees_df$newick[i]
  line2_newick <- trees_df$newick[i + 1]

  tree1 <- ape::read.tree(text = line1_newick)
  tree2 <- ape::read.tree(text = line2_newick)

  spr_dist <- TreeDist::SPRDist(tree1, tree2, method = "deOliveira", symmetric = TRUE)
  info_dist <- TreeDist::MatchingSplitInfoDistance(tree1, tree2)
  rf_dist <- TreeDist::RobinsonFoulds(tree1, tree2)

  results <- rbind(
    results,
    data.frame(start = trees_df$start[i + 1], end = trees_df$end[i+1], spr_dist = spr_dist, splitInfoDist = info_dist, rf_dist = rf_dist)
  )
}

# -----------------------------
# Write output
# -----------------------------
write.table(results, file = output_file, sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)

cat("SPR distances written to:", output_file, "\n")
