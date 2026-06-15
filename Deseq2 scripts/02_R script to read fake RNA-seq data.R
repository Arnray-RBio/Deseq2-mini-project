BiocManager::install("DESeq2")
library(DESeq2)
library(tidyverse)

counts <- read.csv("data/fake_counts.csv", row.names = 1, check.names = FALSE)
metadata <- read.csv("data/fake_metadata.csv")

rownames(metadata) <- metadata$sample_id

counts <- as.matrix(counts)

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = metadata,
  design = ~ condition
)

dds <- dds[rowSums(counts(dds)) >= 10, ]

dds <- DESeq(dds)

res <- results(dds, contrast = c("condition", "treatment", "control"))

res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)

res_df <- res_df %>%
  arrange(padj)

dir.create("results", showWarnings = FALSE)

write.csv(res_df, "results/DESeq2_results_treatment_vs_control.csv", row.names = FALSE)

saveRDS(dds, "results/dds_object.rds")