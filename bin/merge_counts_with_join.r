#!/usr/bin/env Rscript
### merge_counts_with_join.r --- merge counts tables with a full join
##
## Filename: merge_counts_with_join.r
## Author: Zachary Maas <zama8258@colorado.edu>
## Created: Mon Sep 27 10:17:30 2021 (-0600)
##
######################################################################
##
### Commentary:
##
## This script contains code to generate a merged counts from separate
## distinct featureCounts count tables. The primary use cases for this
## are for running featureCounts on individual samples in parallel and
## merging them afterwards to leverage massive parallelism, or to
## combine distinct counts tables from nascent and RNA-seq datasets
## (e.g. without and with exons counted respectively) into one merged
## counts table for things like half life decay analysis. Users should
## note that this script will not preserve both lengths in the latter
## case, using only the lengths and other metadata of the first sample
## provided as input in the command line arguments. This script may
## also be useful for combining counts tables from distinct
## experiments, but that should be done with caution as users should
## ensure that the counts tables were all generated with the same
## parameters.
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

## Load libraries
suppressMessages(library("tidyverse"))
suppressMessages(library("argparse"))

## Build our argument parser
parser <- ArgumentParser()
parser$add_argument("-i", "--input_counts", action="store", dest="input_files",
                    nargs="+", required=TRUE,
                    help="The input counts tables to merge")
parser$add_argument("-o", "--output", action="store", dest="out_file",
                    required=TRUE,
                    help="The file to write the merged counts table to")

## Process arguments
args <- parser$parse_args()
in_files <- args$input_files
out_file <- args$out_file

## Separate data as needed
first_file <- in_files[[1]]
rest_files <- in_files[-1]

## Load and do the merging
merged_dat <- read_delim(first_file, delim='\t')
for (new_file in rest_files) {
    new_dat <- read_delim(new_file, delim='\t') %>%
        select(-c('Chr', 'Start', 'End', 'Strand', 'Length'))
    merged_dat <- full_join(merged_dat, new_dat, by='Geneid')
}

## Write to output file
write_delim(merged_dat, out_file,
            delim='\t', quote='none', col_names=TRUE)

######################################################################
### merge_counts_with_join.r ends here
