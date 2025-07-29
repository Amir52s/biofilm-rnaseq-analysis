# Project: Differential Gene Expression Analysis of P. aeruginosa Biofilm
# Author: AmirHossein Soleymanian
# Date: July, 2025
#
# Description:
# This script analyzes RNA-Seq data from GEO dataset GSE153546 to identify
# genes that are differentially expressed between planktonic and biofilm
# states of Pseudomonas aeruginosa.




# Libraries

library(DESeq2)
library(ggplot2)
library(gprofiler2)


# Loading Data

counts_data <- read.delim("GSE153546_Raw_gene_counts_matrix.txt.gz", header = TRUE, row.names = 1)

# Define the samples for comparison

target_samples <- c("CP1", "CP2", "CP3", "CP4", "CB1", "CB2", "CB3", "CB4")

# Subset the main count matrix to keep only our target samples

counts_subset <- counts_data[, target_samples]

# Creating metadata table (colData)

sample_info <- data.frame(condition = c("Planktonic", "Planktonic", "Planktonic", "Planktonic","Biofilm", "Biofilm", "Biofilm", "Biofilm"), row.names = target_samples)

# DESeq2 Differential Expression Analysis
# First, creating DataSet object

dds <- DESeqDataSetFromMatrix(countData = counts_subset, colData = sample_info, design = ~ condition)

# Run the main DESeq2 analysis

dds <- DESeq(dds)

# Extract the results table, comparing Biofilm vs. Planktonic

res <- results(dds, contrast=c("condition", "Biofilm", "Planktonic"))

# Sort the results by adjusted p-value to see the most significant genes

res_sorted <- res[order(res$padj), ]


# Visualization
# Create an MA-plot for overview

plotMA(res, ylim=c(-3,3), main="MA Plot: Biofilm vs. Planktonic")

# Create a Volcano Plot using ggplot2

res_df <- as.data.frame(res_sorted)
res_df$significant <- ifelse(res_df$padj < 0.05 & !is.na(res_df$padj), "Significant", "Not Significant")

ggplot(data=res_df, aes(x=log2FoldChange, y=-log10(padj), color=significant)) + geom_point(alpha=0.6, size=1.5) + theme_minimal() + scale_color_manual(values=c("Not Significant"="grey50", "Significant"="firebrick")) + labs(title="Volcano Plot: Biofilm vs. Planktonic", x="Log2 Fold Change", y="-Log10 Adjusted P-value") + geom_vline(xintercept=c(-1, 1), linetype="dashed", color="dodgerblue") + geom_hline(yintercept=-log10(0.05), linetype="dashed", color="dodgerblue")


# Final Conclusion

# Get lists of significant genes using a standard adjusted p-value cutoff.

up_genes <- rownames(res_sorted[which(res_sorted$padj < 0.05 & res_sorted$log2FoldChange > 0), ])
down_genes <- rownames(res_sorted[which(res_sorted$padj < 0.05 & res_sorted$log2FoldChange < 0), ])

# Check the number of significant genes
# In this experiment, the numbers were very low (6 up, 4 down).

print(paste("Significant up-regulated genes:", length(up_genes)))
print(paste("Significant down-regulated genes:", length(down_genes)))

# Due to the very low number of significant differentially expressed genes,
# a statistically meaningful enrichment analysis was not possible. 
# The code below is included to show the attempted step.
#
# up_pathways <- gost(query = up_genes, organism = "paeruginosa", sources = "KEGG")
# down_pathways <- gost(query = down_genes, organism = "paeruginosa", sources = "KEGG")
# Final Conclusion:
# The analysis revealed a subtle transcriptional response to biofilm formation
# under these conditions, with only 10 genes being significantly regulated
# (padj < 0.05). This suggests that the transition to a biofilm state in this
# experiment involved specific, targeted changes rather than a large-scale
# reprogramming of major biological pathways.
