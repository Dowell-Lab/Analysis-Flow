#!/usr/bin/env Rscript
### saf_gtf_filter.r --- Filter a GTF file using a SAF file
##
## Filename: saf_gtf_filter.r
## Author: Zachary Maas <zama8258@colorado.edu>
## Created: Mon Feb 15 08:45:09 2021 (-0700)
##
######################################################################
##
### Commentary:
##
## Filter a GTF file using a SAF file. This exists to allow for proper
## exon counting when analyzing RNA-seq data in the pipeline as
## opposed to nascent.
##
######################################################################
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or (at
## your option) any later version.
##
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with GNU Emacs.  If not, see <https://www.gnu.org/licenses/>.
##
######################################################################
##
### Code:

suppressMessages(library("tidyverse"))
suppressMessages(library("argparse"))

## Argument Parsing
parser <- ArgumentParser()
parser$add_argument("-g", "--gtf_file", action="store", dest="gtf_file",
                    help="The gtf table to use")
parser$add_argument("-s", "--saf_file", action="store", dest="saf_file",
                    help="The saf table to use")
parser$add_argument("-o", "--out_file", action="store", dest="out_file",
                    help="The out file to use")
args <- parser$parse_args()

## Sample files (TODO switch to argparse)
gtf <- args$gtf_file
saf <- args$saf_file
out_gtf <- args$out_file

gtf_df <- read_delim(gtf, delim="\t", col_names=FALSE) %>%
    separate(X9, c("n1", "geneid", "n2", "txid"), sep = " ") %>%
    mutate(geneid = str_remove_all(geneid, '[;"]'))

saf_df <- read_delim(saf, delim="\t", col_names=FALSE) %>%
    mutate(geneid = X1)

filtered_df <- gtf_df %>% semi_join(saf_df, by=c("geneid")) %>%
    mutate(geneid=paste0('"', geneid, '";')) %>%
    unite(metadata, n1, geneid, n2, txid, sep=" ")

write_tsv(filtered_df, out_gtf, col_names = FALSE, quote_escape = "none")

######################################################################
### saf_gtf_filter.r ends here
