library(tximport)
library(DESeq2)

samples <- read.csv("samples.csv", stringsAsFactors = FALSE)
samples$condition <- factor(samples$condition, levels = c("control", "mazF20m"))

files <- file.path("quant", samples$sample, "quant.sf")
names(files) <- samples$sample
stopifnot(all(file.exists(files)))

txi <- tximport(files, type = "salmon", txOut = TRUE, dropInfReps = TRUE)

cat("\n--- Raw library sizes (mapped reads) ---\n")
print(round(colSums(txi$counts)))

dds <- DESeqDataSetFromTximport(txi, colData = samples, design = ~ condition)

keep <- rowSums(counts(dds) >= 10) >= 3
cat("\nTranscripts before filtering:", nrow(dds), "\n")
dds <- dds[keep, ]
cat("Transcripts after filtering:", nrow(dds), "\n")

dds <- DESeq(dds)

cat("\n--- Normalization factors (geometric mean per sample) ---\n")
nf <- normalizationFactors(dds)
print(round(exp(colMeans(log(nf))), 3))

res <- results(dds, contrast = c("condition", "mazF20m", "control"))
res <- res[order(res$padj), ]

cat("\n--- Summary ---\n")
summary(res)

cat("\nPadj < 0.05:", sum(res$padj < 0.05, na.rm = TRUE), "\n")
cat("Padj < 0.05 & |LFC| > 1:", sum(res$padj < 0.05 & abs(res$log2FoldChange) > 1, na.rm = TRUE), "\n")

write.csv(as.data.frame(res), "results/deseq2_results.csv")

vsd <- vst(dds, blind = TRUE)
png("results/PCA.png", width = 800, height = 600)
print(plotPCA(vsd, intgroup = "condition"))
dev.off()

cat("\n--- Top 15 by padj ---\n")
print(head(as.data.frame(res)[, c("baseMean","log2FoldChange","pvalue","padj")], 15))
