#!/usr/bin/env python3
"""Script to split lines into 100 for metagene"""

import argparse
import math

## Set up argument parsing
PARSER = argparse.ArgumentParser()
PARSER.add_argument("-f", "--file")
PARSER.add_argument("-n", "--num_regions")
PARSER.add_argument("-o", "--out")
ARGS = PARSER.parse_args()

FILE = ARGS.file
NUM_REGIONS = int(ARGS.num_regions)
OUTFILE = ARGS.out
BED = False

if NUM_REGIONS is None or FILE is None or OUTFILE is None:
    print("Missing arguments. Try --help")
    quit()

with open(OUTFILE, "w+") as out:
    # out.write('\t'.join(["GeneID", "Chr", "Start", "End", "Strand"]) + '\n')
    with open(FILE) as f:
        for line in f:
            curr_line = line.split()
            newline = line.strip()
            spread = math.floor(
                (int(curr_line[2]) - int(curr_line[1])) / NUM_REGIONS
            )
            for i in range(NUM_REGIONS):
                # Featurecounts partial - Fractional Overlap
                c_start = str(int(curr_line[1]) + (i * spread))
                c_end = str(int(curr_line[1]) + ((i + 1) * spread))
                c_id = curr_line[3] + "/" + str(i)
                if BED is True:
                    # BED Format
                    out.write(
                        "\t".join(
                            [
                                curr_line[0],
                                c_start,
                                c_end,
                                c_id,
                                curr_line[4],
                                curr_line[5],
                            ]
                        )
                        + "\n"
                    )
                else:
                    # SAF Format for featureCounts
                    out.write(
                        "\t".join(
                            [c_id, curr_line[0], c_start, c_end, curr_line[5]]
                        )
                        + "\n"
                    )
