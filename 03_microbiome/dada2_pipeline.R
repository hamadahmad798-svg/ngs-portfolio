library(dada2)

path <- "data/MiSeq_SOP"
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names=TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names=TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
cat("Samples found:", length(sample.names), "\n")
print(sample.names)

filtFs <- file.path("results/filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path("results/filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names; names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     truncLen=c(240,160), maxN=0, maxEE=c(2,2),
                     truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=4)
cat("\n--- Filtering ---\n"); print(head(out))

cat("\n--- Learning error rates ---\n")
errF <- learnErrors(filtFs, multithread=4)
errR <- learnErrors(filtRs, multithread=4)

dadaFs <- dada(filtFs, err=errF, multithread=4)
dadaRs <- dada(filtRs, err=errR, multithread=4)

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=FALSE)
seqtab <- makeSequenceTable(mergers)
cat("\nASVs before chimera removal:", ncol(seqtab), "\n")
cat("Amplicon length distribution:\n"); print(table(nchar(getSequences(seqtab))))

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=4, verbose=TRUE)
cat("ASVs after chimera removal:", ncol(seqtab.nochim), "\n")
cat("Fraction of reads retained:", round(sum(seqtab.nochim)/sum(seqtab), 4), "\n")

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN),
               sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input","filtered","denoisedF","denoisedR","merged","nonchim")
rownames(track) <- sample.names
cat("\n--- Read tracking ---\n"); print(track)

saveRDS(seqtab.nochim, "results/seqtab_nochim.rds")
