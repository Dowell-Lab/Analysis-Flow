# Analysis Flow Pipeline

**Pipelines do not supercede scientific critical thinking. You need to
understand exactly what analyses you're doing and why.**

You should not run this pipeline without looking at the code for the
analysis you're interested in and understanding the pitfalls and sharp
edges of the tools used. This will not automatically analyze your data
and then tell you exactly how to think about the results. The analyses
in this pipeline are those commonly used in publications, and
effectively generating those analyses and interpreting the results
requires not just an understanding of the tools used but also the
biological context of your experiment.

This repository contains Zach's Analysis-Flow pipeline to perform
basic analysis on processed sequencing data (Nascent and RNA-Seq).
This pipeline performs the following:
- Isoform filtering over RefSeq Genes
- Differential Expression Analysis
- Pause Index Analysis
- Principal Component Analysis
- Metagene Plot Generation

This pipeline can also differentiate between single and paired end
data as well as Nascent and RNA-seq protocols using a design table.
This allows for complex experimental designs with disparate types of
data to be analyzed together.

To use this pipeline, you need to have run the Dowell Lab's
RNASeq-Flow or Nascent-Flow pipeline first. Then, edit the
configuration file in `conf/example.conf` with the paths of your
processed data and any parameters that are organism dependent.

You will also need to install the following binary dependencies:
``` shell
subread
bedtools
```
And the following R depndencies:
``` shell
tidyverse
ggplot2
ggfortify
plyr
argparse
ggthemes
reshape2
digest
limma
sva
DESeq2
```

You can easily install the R dependencies with the following commands:
``` shell
R -e "install.packages(c('tidyverse', 'ggplot2', 'ggfortify', 'plyr', 'argparse', 'ggthemes', 'reshape2', 'digest'), repos='http://cran.rstudio.com/')"
R -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')"
R -e "BiocManager::install(c('limma', 'sva', 'DESeq2'))"
```

## Experimental Features

This pipeline currently supports some experimental features for making
both nascent and rna-seq data run successfully in a single run. Ask
Zach for more detail.

## Incomplete Tasks

- Wrap up all the dependencies in a Singularity container to make
  running the pipeline a zero-effort kind of thing.
- Add a final step that does automatic report generation for
  end-users. This involves getting all their figures into a single zip
  archive along with a document on how to interpret the results.
- Clean up inputs to reduce the amount of file copying and speed up R scripts.

## License

Analysis-Flow Pipeline
Copyright (C) 2020 Dowell Lab

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
