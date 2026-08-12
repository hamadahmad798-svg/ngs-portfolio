library(tximport); library(DESeq2)
samples <- read.csv("samples.csv", stringsAsFactors = FALSE)
samples$condition <- factor(samples$condition, levels = c("control","mazF20m"))
files <- file.path("quant", samples$sample, "quant.sf")
names(files) <- samples$sample
txi <- tximport(files, type="salmon", txOut=TRUE, dropInfReps=TRUE)
dds <- DESeqDataSetFromTximport(txi, colData=samples, design=~condition)
dds <- dds[rowSums(counts(dds) >= 10) >= 3, ]
vsd <- vst(dds, blind=TRUE)
pc <- plotPCA(vsd, intgroup="condition", returnData=TRUE)
pc$libsize <- colSums(txi$counts)[rownames(pc)]
print(pc[, c("PC1","PC2","condition","libsize")])
