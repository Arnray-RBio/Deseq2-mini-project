library(DESeq2)
library(tidyverse)
library(pheatmap)
library(ggrepel)

# Read saved DESeq2 object and results
dds <- readRDS("results/dds_object.rds")

metadata <- read.csv("data/fake_metadata.csv")
rownames(metadata) <- metadata$sample_id

res_df <- read.csv("results/DESeq2_results_treatment_vs_control.csv")

# Create figures folder
dir.create("figures", showWarnings = FALSE)

# -----------------------------
# 1. PCA plot
# -----------------------------

vsd <- vst(dds, blind = FALSE)

pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw() +
  labs(title = "PCA plot")

ggsave("figures/PCA_plot.pdf", p_pca, width = 6, height = 5)


# -----------------------------
# 2. Volcano plot
# -----------------------------

volcano_df <- res_df %>%
  filter(!is.na(padj)) %>%
  mutate(
    significance = case_when(
      padj < 0.05 & log2FoldChange > 1 ~ "Upregulated",
      padj < 0.05 & log2FoldChange < -1 ~ "Downregulated",
      TRUE ~ "Not significant"
    )
  )

p_volcano <- ggplot(volcano_df, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
  geom_point(alpha = 0.7) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_bw() +
  labs(
    x = "log2 fold change",
    y = "-log10 adjusted p-value",
    title = "Volcano plot: treatment vs control"
  )

ggsave("figures/volcano_plot.pdf", p_volcano, width = 7, height = 5)


# -----------------------------
# 3. Heatmap of top genes
# -----------------------------

top_genes <- res_df %>%
  filter(!is.na(padj)) %>%
  arrange(padj) %>%
  dplyr::slice_head(n = 30) %>%
  pull(gene_id)

mat <- assay(vsd)[top_genes, ]

# Scale each gene across samples
mat_scaled <- t(scale(t(mat)))

annotation_col <- metadata["condition"]

pdf("figures/top_gene_heatmap.pdf", width = 7, height = 8)

pheatmap(
  mat_scaled,
  annotation_col = annotation_col,
  show_rownames = TRUE,
  show_colnames = TRUE,
  main = "Top differentially expressed genes"
)

dev.off()