library(tidyverse)

set.seed(123)

metadata <- data.frame(
  sample_id = paste0("sample_", 1:6),
  condition = rep(c("control", "treatment"), each = 3)
)

rownames(metadata) <- metadata$sample_id

gene_ids <- paste0("gene_", 1:1000)

counts <- matrix(
  rnbinom(1000 * 6, mu = 100, size = 1),
  nrow = 1000,
  ncol = 6
)

rownames(counts) <- gene_ids
colnames(counts) <- metadata$sample_id

counts[1:25, metadata$condition == "treatment"] <- counts[1:25, metadata$condition == "treatment"] * 4
counts[26:50, metadata$condition == "treatment"] <- round(counts[26:50, metadata$condition == "treatment"] * 0.25)

dir.create("data", showWarnings = FALSE)

write.csv(counts, "data/fake_counts.csv")
write.csv(metadata, "data/fake_metadata.csv", row.names = FALSE)

