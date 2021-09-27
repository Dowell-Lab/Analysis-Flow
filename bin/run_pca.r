#!/usr/bin/env Rscript
### run_pca_r.r --- run pca with batch correction
##
## Filename: run_pca_r.r
## Description: Run PCA
## Author: Student Zachary Maas <zama8258@colorado.edu>
## Maintainer: Student Zachary Maas <zama8258@colorado.edu>
## Created: Thu Jan 31 14:25:39 2019 (-0700)
##
######################################################################
##
### Commentary:
##
## Runs PCA with batch correction...
##
### Code:

suppressMessages(library("tidyverse"))
suppressMessages(library("ggplot2"))
suppressMessages(library("ggfortify"))
suppressMessages(library("plyr"))
suppressMessages(library("argparse"))
suppressMessages(library("limma"))
suppressMessages(library("sva"))

parser <- ArgumentParser()
parser$add_argument("-c", "--counts", action="store", dest="countsFile",
                    help="The counts table to use for PCA")
args <- parser$parse_args()

countsFile <- args$countsFile

counts <- read_delim(countsFile, delim = "\t")

tpm_normalize <- function(data, sample) {
    ## First, calculate RPK for the strand only
    RPK_sample <- (data[[sample]] / (data$Length / (10 ^ 3)))
    ## Then, calculate RPK for both strands (for scaling factor)
    RPK <- ((data[[sample]]) /
            (data$Length / (10 ^ 3)))
    ## Then, calculate the scaling factor
    scale <- sum(RPK) / 1000000
    ## Divide RPK values by scaling factor
    out <- RPK_sample / scale
    return(out)
}

counts_orig <- counts
counts <- counts %>%
    subset(select = -c(Geneid, Chr, Start, End, Strand, Length))

samples <- colnames(counts)
fx <- data.frame(lapply(samples, function(x) tpm_normalize(counts_orig, x)))
colnames(fx) <- samples

counts_t <- data.frame(t(fx))
pca_raw <- prcomp(counts_t)
pcar <- data.frame(pca_raw$x)
pcar$group <- samples

ggplot() +
    geom_point(data = pcar, mapping = aes(x = PC1, y = PC2,
                                          color = samples, size = 6)) +
    labs(x = "PC1", y = "PC2",
         title = "Principal Component Analysis of Samples",
         color = "Condition") +
    guides(size = FALSE)
ggsave(str_c("pca.pdf"),
       plot = last_plot(),
       height = 5,
       width = 10)

######################################################################
### run_pca_r.r ends here
