#!/usr/bin/env Rscript
### metagene_graph_custom.r --- Custom Metagene Plots
##
## Filename: metagene_graph_custom.r
## Author: Student Zachary Maas <zama8258@colorado.edu>
## Maintainer: Student Zachary Maas <zama8258@colorado.edu>
##
######################################################################
##
### Commentary:
##
## Builds Metagene Plots from Custom CountsSense File
##
######################################################################
##
### Code:

suppressMessages(library("tidyverse"))
suppressMessages(library("ggplot2"))
suppressMessages(library("argparse"))
suppressMessages(library("reshape2"))
suppressMessages(library("digest"))

parser <- ArgumentParser()
parser$add_argument("-s", "--sense", action="store", dest="countsSense",
                    help="The sense reads counts table to use")
parser$add_argument("-a", "--antisense", action="store", dest="countsAntiSense",
                    help="The antisense reads table to use")
parser$add_argument("-gi", "--group_i", action="store", dest="group_i",
                    nargs="+",
                    help="The members of the first group.")
parser$add_argument("-gj", "--group_j", action="store", dest="group_j",
                    nargs="+",
                    help="The members of the second group.")
parser$add_argument("-ni", "--name_i", action="store", dest="name_i",
                    help="The name of the first group.")
parser$add_argument("-nj", "--name_j", action="store", dest="name_j",
                    help="The name of the second group.")
parser$add_argument("-n", "--numbins", action="store", dest="numbins",
                    help="The number of bins used")
parser$add_argument("-o", "--output", action="store", dest="outfile",
                    help="The output image file name.")
args <- parser$parse_args()

sense <- args$countsSense
antisense <- args$countsAntiSense
group_i <- args$group_i
group_j <- args$group_j
name_i <- args$name_i
name_j <- args$name_j
num_bins <- as.numeric(args$numbins)
outfile <- args$outfile

countsSense <- read_delim(sense, delim="\t")
countsAntiSense <- read_delim(antisense, delim="\t")

##setwd("/home/zach/dowell_lab/pausing_meta_analysis/data")
##countsSense <- read_delim('metagene_counts_sense_fix.txt', delim="\t")
##countsAntiSense <- read_delim('metagene_counts_antisense_fix.txt', delim="\t")
##num_bins <- 100

df_sense <- countsSense %>% separate(Geneid, into = c("geneid", "coord"), sep="/")
df_sense$coord = as.numeric(df_sense$coord)
df_antisense <- countsAntiSense %>% separate(Geneid, into = c("geneid", "coord"), sep = "/")
df_antisense$coord = as.numeric(df_antisense$coord)

print("Using Groups:")
print(name_i)
group_i <- lapply(group_i, function(x) paste0(x, ".sorted.bam"))
print(group_i)
print(name_j)
group_j <- lapply(group_i, function(x) paste0(x, ".sorted.bam"))
print(group_j)

meltify_frame <- function(frame) {
    new_frame <- frame %>%
        subset(select = -c(geneid, Chr, Start, End)) %>%
        melt(id = c("coord", "Strand", "Length")) %>%
        mutate(condition = ifelse(is.element(variable, group_i), name_i, name_j),
               coord = ifelse(Strand == '-', (num_bins - 1) - coord, coord),
               reads = value) %>%
        subset(select = -c(variable, value)) %>%
        as_tibble()
    return(new_frame)
}

sense_melt <- meltify_frame(df_sense)
antisense_melt <- meltify_frame(df_antisense)

## Normalize the counts for each region by TPM
tpm_normalize <- function(fst, snd, sample) {
    ## First, calculate RPK for the strand only
    RPK_sample <- (fst[[sample]] / (fst$Length / (10 ^ 3)))
    ## Then, calculate RPK for both strands (for scaling factor)
    RPK <- ((fst[[sample]] + snd[[sample]]) /
            (fst$Length / (10 ^ 3)))
    ## Then, calculate the scaling factor
    scale <- sum(RPK) / 1000000
    ## Divide RPK values by scaling factor
    out <- RPK_sample / scale
    return(out)
}

sense_melt$reads <- tpm_normalize(sense_melt, antisense_melt, "reads")
antisense_melt$reads <- -tpm_normalize(antisense_melt, sense_melt, "reads")

summary_stats <- function(frame) {
    new_frame <- frame %>% group_by(coord, condition) %>%
        summarize(mu = mean(reads),
                  var = sd(reads) / sqrt(length(reads))) %>%
        mutate(min = mu - var, max = mu + var)
    return(new_frame)
}

sense_melt_final <- summary_stats(sense_melt)
antisense_melt_final <- summary_stats(antisense_melt)

hex_to_int <- function(h) {
    xx = strsplit(tolower(h), "")[[1L]]
    pos = match(xx, c(0L:9L, letters[1L:6L]))
    sum((pos - 1L) * 16^(rev(seq_along(xx) - 1)))
}
sample_color <- function(name) {
    hue <- (hex_to_int(digest(name, algo='xxhash32')) %% 1009) / 1009
    new_color <- hsv(hue, 1, 0.75)
    return(new_color)
}

library('ggthemes')
c(sample_color(name_i), sample_color(name_j))
ggplot() + theme_tufte() +
    scale_color_manual(values =  c(sample_color(name_i), sample_color(name_j))) +
    scale_fill_manual(values =  c(sample_color(name_i), sample_color(name_j))) +
    geom_line(data = sense_melt_final, aes(x = coord, y = mu, color = condition)) +
    geom_line(data = antisense_melt_final, aes(x = coord, y = mu, color = condition)) +
    geom_ribbon(data = sense_melt_final, aes(x = coord,
                                             ymin = min, ymax = max,
                                             fill = condition), alpha = 0.2) +
    geom_ribbon(data = antisense_melt_final, aes(x = coord,
                                                 ymin = min, ymax = max,
                                                 fill= condition), alpha = 0.2) +
    geom_hline(yintercept = 0) +
    labs(title = paste0(name_i, " vs ", name_j),
         x = "Relative Position", y = "Read Depth",
         color = "Condition", fill = "SD of Mean") +
    ggsave(outfile, width = 16, height = 9)
## Output Data Frame DataFrame

######################################################################
### metagene_graph_custom.r ends here
