seqtab <- readRDS("results/seqtab_nochim.rds")
taxa   <- readRDS("results/taxa.rds")

bad <- (taxa[,"Family"] %in% "Mitochondria") | (taxa[,"Order"] %in% "Chloroplast")
cat("ASVs flagged as host/plastid:", sum(bad, na.rm=TRUE), "\n")
cat("Reads removed:", sum(seqtab[, which(bad)]), "of", sum(seqtab),
    sprintf("(%.3f%%)\n", 100*sum(seqtab[, which(bad)])/sum(seqtab)))

seqtab <- seqtab[, !bad]
taxa   <- taxa[!bad, ]
saveRDS(seqtab, "results/seqtab_clean.rds")
saveRDS(taxa,   "results/taxa_clean.rds")
cat("ASVs retained:", ncol(seqtab), "\n")
