#! /usr/bin/env Rscript

library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
    stop("Usage: Rscript plot_depth.R <depth.tsv> <tophit.tsv> <output.png>")
}

depth_file <- args[1]
tophit_file <- args[2]
output_file <- args[3]

depth <- read.table(depth_file, header = FALSE, comment.char = "",
                    col.names = c("reference", "position", "depth"))

tophit <- read.table(tophit_file, header = FALSE, comment.char = '',
                    col.names = c("sample","segment","reference","mapped","coverage","meandepth"))

order <- tophit$reference

depth <- subset(depth, reference %in% order)

depth$reference <- factor(depth$reference, levels = order)

p <- ggplot(depth, aes(position, depth)) +
    geom_line() +
    facet_wrap(~reference, scales = "free_x", ncol = 1) +
    theme_bw() +
    labs(x = "Position", y = "Depth")

ggsave(
    filename = output_file,
    plot = p,
    width = 8,
    height = 11,
    dpi = 300
)