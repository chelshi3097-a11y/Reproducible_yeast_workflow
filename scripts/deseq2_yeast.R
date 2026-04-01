# ================================
# RNA-seq DESeq2 Analysis - Yeast
# ================================

set.seed(42)

library(DESeq2)
library(ggplot2)
library(pheatmap)
library(ggrepel)

# -------------------------------
# Step 1: File paths from command line
# -------------------------------

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 8) {
  stop("Expected 8 arguments: counts metadata results ma volcano pca heatmap sampledist")
}


counts_file   <- args[1]
metadata_file <- args[2]
results_file  <- args[3]
ma_file       <- args[4]
volcano_file  <- args[5]
pca_file      <- args[6]
heatmap_file  <- args[7]
sampledist_file <- args[8]



# -------------------------------
# Step 2: Read featureCounts
# -------------------------------
fc <- read.table(
  counts_file,
  header = TRUE,
  sep = "\t",
  comment.char = "#",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

count_data <- fc[, 7:ncol(fc)]
rownames(count_data) <- fc$Geneid

clean_names <- function(x) {
  x <- gsub(".*/", "", x) # remove folder path
  x <- gsub("\\.Aligned.sortedByCoord\\.out\\.bam$", "", x) # remove suffix correctly
  x
}

colnames(count_data) <- clean_names(colnames(count_data))

# -------------------------------
# Step 3: Read metadata
# -------------------------------
coldata <- read.csv(metadata_file, stringsAsFactors=FALSE)

# Use Run_Accession as rownames (must match count matrix colnames)
rownames(coldata) <- coldata$Run_Accession
coldata$condition <- factor(coldata$Group, levels = c("WT", "MUT"))

# -------------------------------
# Step 4: Match samples safely
# -------------------------------
missing <- setdiff(rownames(coldata), colnames(count_data))
if (length(missing) > 0) {
  stop("Samples missing in count matrix: ", paste(missing, collapse = ", "))
}

count_data <- count_data[, rownames(coldata)]
stopifnot(identical(rownames(coldata), colnames(count_data)))

# -------------------------------
# Step 5: DESeq2 object
# -------------------------------
dds <- DESeqDataSetFromMatrix(
  countData = count_data,
  colData   = coldata,
  design    = ~ condition
)

dds <- DESeq(dds)
res <- results(dds)

res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)

res_df <- res_df[!is.na(res_df$padj), ]

# avoid -log10(0)
if (any(res_df$padj == 0)) {
  res_df$padj[res_df$padj == 0] <- min(res_df$padj[res_df$padj > 0])
}

write.csv(res_df, results_file, row.names = FALSE)


# -------------------------------
# Step 6: MA plot
# -------------------------------
png(ma_file)

plotMA(res, ylim = c(-5, 5), main = "MA Plot")
dev.off()

# -------------------------------
# Step 7: Volcano plot
# -------------------------------
res_df$significant <- ifelse(res_df$padj < 0.05, "yes", "no")
sig_genes <- res_df[res_df$significant == "yes", ]

genes_to_label <- data.frame()
if (nrow(sig_genes) > 0) {
  top_up   <- head(sig_genes[order(-sig_genes$log2FoldChange), ], 5)
  top_down <- head(sig_genes[order(sig_genes$log2FoldChange), ], 5)
  genes_to_label <- rbind(top_up, top_down)
}

volcano_plot <- ggplot(
  res_df,
  aes(x = log2FoldChange, y = -log10(padj), color = significant)
) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("no" = "grey", "yes" = "red")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  theme_minimal(base_size = 12) +
  labs(title = "Volcano Plot (padj < 0.05)")

if (nrow(genes_to_label) > 0) {
  volcano_plot <- volcano_plot +
    geom_text_repel(
      data = genes_to_label,
      aes(label = gene),
      size = 3
    )
}

ggsave(volcano_file, volcano_plot, width=7, height=5, dpi=300)


# -------------------------------
# Step 8: VST + PCA (FIXED)
# -------------------------------
vsd <- vst(dds, blind = FALSE)

pca_df <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- attr(pca_df, "percentVar")

pca_plot <- ggplot(
  pca_df,
  aes(x = PC1, y = PC2, color = condition, label = name)
) +
  geom_point(size = 3) +
  geom_text_repel(size = 4) +
  xlab(paste0("PC1: ", round(percentVar[1] * 100, 1), "%")) +
  ylab(paste0("PC2: ", round(percentVar[2] * 100, 1), "%")) +
  theme_minimal(base_size = 14) +
  ggtitle("PCA Plot of Samples")

ggsave(pca_file, pca_plot, width=7, height=5, dpi=300)


# -------------------------------
# Step 9: Top 20 DE genes heatmap
# -------------------------------
n_top <- min(20, nrow(res_df))
top_genes <- res_df[order(res_df$padj), ][1:n_top, ]

top_mat <- assay(vsd)[top_genes$gene, , drop = FALSE]
top_mat_scaled <- t(scale(t(top_mat)))

png(heatmap_file, 800, 600)

pheatmap(
  top_mat_scaled,
  annotation_col = coldata["condition"],
  fontsize_row = 8,
  main = "Top DE Genes"
)
dev.off()

# -------------------------------
# Step 10: Sample distance heatmap
# -------------------------------
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)

png(sampledist_file, 800, 600)

pheatmap(
  sampleDistMatrix,
  annotation_col = coldata["condition"],
  main = "Sample-to-Sample Distance"
)
dev.off()
