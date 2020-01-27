## Script to run DESeq2
## Maintainer: Zachary Maas <zama8258@colorado.edu>

#######################
## Preliminary Setup ##
#######################

suppressMessages(library("tidyverse"))
suppressMessages(library("ggthemes"))
suppressMessages(library("DESeq2"))
suppressMessages(library("argparse"))

## Argument Parsing
parser <- ArgumentParser()
parser$add_argument("-c", "--counts_file", action="store", dest="counts_file",
                    help="The counts table to use")
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
parser$add_argument("-t", "--conversion_table", action="store", dest="conversion_file",
                    help="The gene id conversion to use.")
parser$add_argument("-o", "--outdir", action="store", dest="outdir",
                    help="The output directory name.")
args <- parser$parse_args()

counts_file <- args$counts_file
group_i <- args$group_i
group_j <- args$group_j
name_i <- args$name_i
name_j <- args$name_j
conversion_file <- args$conversion_file
outdir <- args$outdir
comparison <- paste0(name_i, "_vs_", name_j)
legend_title <- paste0("(", name_i, " / ", name_j, ")")

setwd(outdir)
print("parsed args / imported data")

## Set up a function for plotting
## This is used later to do all the analysis + plotting
makefig <- function(deseqdata, fileprefix) {
    res <- DESeq2::results(deseqdata)
    table(res$padj < 0.05)
    ## Order by adjusted p-value
    res <- res[order(res$padj), ]
    ## Merge with normalized count data
    resdata <- merge(as.data.frame(res), as.data.frame(counts(deseqdata, normalized = TRUE)),
                     by = "row.names", sort = FALSE)
    names(resdata)[1] <- "Gene"
    head(resdata)
    ## Write results in the default format
    write.csv(resdata, file = paste0(fileprefix, "-diffexpr-results.csv"))

                                        # Generate the MA plot
    p_thresh <- 0.005
    p_str_signif <- str_c("p < ", p_thresh)
    p_str_not_signif <- str_c("p >= ", p_thresh)
    black <- "#1b1b1b"
    red <- "#d50000"
    res$pvalue <- res$pvalue %>% replace_na(1)
    res_df <- as_tibble(res) %>% mutate(signif = pvalue < p_thresh)
    ylims <- c(quantile(res_df$log2FoldChange, 0.001, na.rm = TRUE),
               quantile(res_df$log2FoldChange, 0.999, na.rm = TRUE))
    ggplot(data = res_df) +
        geom_point(aes(x = baseMean, y = log2FoldChange, color = signif),
                   size = 1.00) +
        scale_x_log10() +
        lims(y = ylims) +
        scale_color_manual(values = c(black, red),
                           name="Significance",
                           labels=c(p_str_not_signif,
                                    p_str_signif)) +
        theme_tufte() +
        geom_hline(yintercept = 0, color = red) +
        theme(legend.position=c(0.9,0.9)) +
        theme(legend.background = element_rect(colour = "#eceff1",
                                               fill = "white", linetype='solid')) +
        labs(title = paste0("DESeq2 Differential Expression ", legend_title),
             x = "Mean of Normalized Counts",
             y = paste0("Log2 Fold-Change ", legend_title)) +
        ggsave(paste0(fileprefix, "-diffexpr-maplot.pdf"), width = 8, height = 4.5)

    ## Generate a tsv with results
    write.table(na.omit(as.data.frame(res[2])),
                file = paste0(fileprefix, "-genes_for_gsea.txt"), row.names = TRUE,
                quote = FALSE)
    write.table(res, file = paste0(fileprefix, "-DESeq.res.txt"), append = FALSE, sep = "\t")

    ## Generate preranked file for GSEA
    rnkdf <- tibble(gene = rownames(res), rnk = res$pvalue / sign(res$log2FoldChange)) %>%
        arrange(desc(rnk)) %>% drop_na()
    write.table(rnkdf, file = paste0(fileprefix, ".rnk"),
                append = FALSE, col.names = FALSE, row.names = FALSE,
                quote = FALSE, sep = "\t")

    return(res)
}

## ID Conversion Table
idConvert <- read_delim(conversion_file,
                        col_names = c("rowname", "common"), delim = " ")

## Experimental condition
condition <- factor(c(rep(name_i, length(group_i)), rep(name_j, length(group_j))),
                    levels=c(name_i, name_j))
condition <- relevel(condition, name_j)
print(condition)

## We have to fix the counts table by removing the first row, hence "counts_fix.txt"
seqdata <- read.delim(counts_file, stringsAsFactors = FALSE, header = TRUE,
                           row.names = 1)
## Then, we filter to only include the resolved isoforms
df <- as_tibble(rownames_to_column(seqdata)) %>%
    left_join(idConvert) %>%
    distinct(common, .keep_all = TRUE)

## Set the rownames like DESeq2 requires
dt <- data.frame(df)
rownames(dt) <- dt$common
dt$common <- NULL

## Now, actually perform DE-Seq
countdata <- dt[, 7:ncol(dt)]
countdata <- as.matrix(countdata)

coldata <- data.frame(row.names = colnames(countdata), condition)
countdata <- countdata[, rownames(coldata)]

if (!all(rownames(coldata) == colnames(countdata))) {
    print("Error. Matrix and count labels do not match.")
    quit(status = 1, "no")
}

dds <- DESeqDataSetFromMatrix(countData = countdata,
                              colData = coldata,
                              design = ~condition)
dds <- DESeq(dds)

ctrl <- makefig(dds, comparison)
