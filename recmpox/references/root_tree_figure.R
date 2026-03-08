#!/usr/bin/env Rscript
#
# Midpoint-root a phylogeny and plot tip labels as circles with clade colors.
# Writes rooted.tree and tree_figure.pdf in the same directory as the input tree.
#
# Usage: Rscript root_tree_figure.R [treefile]
# Default: alignment.treefile (when run from phylogeny dir)
#
# Requires: ape, phytools, ggtree, ggplot2.

args <- commandArgs(trailingOnly = TRUE)
tree_path <- if (length(args) >= 1) args[1] else "alignment.treefile"
out_dir <- dirname(tree_path)
if (out_dir == "") out_dir <- "."

if (!file.exists(tree_path)) {
  stop("Tree file not found: ", tree_path)
}

suppressPackageStartupMessages({
  library(ape)
  library(phytools)
  library(ggtree)
  library(ggplot2)
})
old_warn <- getOption("warn")
on.exit(options(warn = old_warn), add = TRUE)
options(warn = -1)

tr <- read.tree(tree_path)
tr <- phytools::midpoint.root(tr)
rooted_path <- file.path(out_dir, "rooted.tree")
write.tree(tr, rooted_path)

labels <- tr$tip.label
type <- ifelse(grepl("_tracts", labels, fixed = TRUE), "recombinant", "reference")
clade <- rep("Ia", length(labels))
clade[type == "recombinant"] <- "Recombinant"
clade[grepl("sh2024Ib", labels, fixed = TRUE)] <- "sh2024Ib"
clade[grepl("sh2023Ib", labels, fixed = TRUE)] <- "sh2024Ib"
clade[grepl("sh2024[iI]a", labels)] <- "sh2024Ia"
clade[grepl("sh2017IIb", labels, fixed = TRUE)] <- "sh2017IIb"
clade[grepl("_IIa_", labels, fixed = TRUE)] <- "IIa"
clade[grepl("^IIa_", labels)] <- "IIa"
tip_data <- data.frame(label = labels, type = type, clade = clade, stringsAsFactors = FALSE)
tip_data$plot_clade <- ifelse(tip_data$clade == "sh2024Ib", "sh2023Ib", tip_data$clade)

p <- ggtree(tr, linewidth = 0.65)
clade_colors <- c(
  Recombinant = "#FF6B9D",
  Ia = "#1e5f72",
  sh2024Ia = "#b83c28",
  sh2023Ib = "#2d7a4a",
  sh2017IIb = "#c97a08",
  IIa = "#5c2270"
)
tip_size <- 4.2
legend_order <- c("sh2023Ib", "sh2024Ia", "Ia", "sh2017IIb", "IIa", "Recombinant")
legend_labels <- c(
  sh2023Ib = "Ib (sh2023Ib)",
  sh2024Ia = "Ia (sh2024Ia)",
  Ia = "Ia (non-sh2024Ia)",
  sh2017IIb = "IIb (sh2017IIb)",
  IIa = "IIa",
  Recombinant = "Potential recombinant"
)

p <- p %<+% tip_data +
  geom_tippoint(aes(fill = plot_clade), size = tip_size, shape = 21, color = "black", stroke = 0.5) +
  scale_fill_manual(
    values = clade_colors,
    breaks = legend_order,
    labels = legend_labels,
    na.value = "grey70",
    name = "Clade",
    drop = FALSE
  ) +
  theme(
    legend.position = c(0.18, 0.97),
    legend.justification = c(0.5, 1),
    legend.background = element_rect(fill = "white", colour = "black", linewidth = 0.4),
    legend.margin = margin(5, 5, 5, 5),
    legend.spacing.y = unit(4.5, "mm"),
    legend.key.size = unit(6.75, "mm"),
    legend.title = element_text(face = "bold", size = 15),
    legend.text = element_text(size = 13.5),
    text = element_text(family = "sans")
  )

tree_depth <- max(ape::node.depth.edgelength(tr))
p <- p +
  geom_tiplab(aes(label = label), size = 2.4, hjust = -0.04, align = TRUE, linesize = 0.2) +
  xlim_tree(tree_depth * 1.25)
p <- p + geom_treescale(x = 0, y = -1, width = 0.002, fontsize = 2.8, linesize = 0.4)

fig_height <- max(8, length(tr$tip.label) * 0.15)
pdf_path <- file.path(out_dir, "tree_figure.pdf")
ggsave(pdf_path, p, width = 12, height = fig_height, limitsize = FALSE)
# Use base R SVG device (no svglite package required)
svg_path <- file.path(out_dir, "tree_figure.svg")
ggsave(svg_path, p, width = 12, height = fig_height, limitsize = FALSE, device = grDevices::svg)
