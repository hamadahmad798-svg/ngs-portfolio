seqtab <- readRDS("results/seqtab_nochim.rds")
taxa   <- readRDS("results/taxa.rds")

cat("--- Assignment success by rank ---\n")
n <- nrow(taxa)
for (r in colnames(taxa)) {
  ok <- sum(!is.na(taxa[, r]))
  cat(sprintf("%-8s %4d / %4d  (%5.1f%%)\n", r, ok, n, 100*ok/n))
}

cat("\n--- Phylum distribution (ASV counts) ---\n")
print(sort(table(taxa[,"Phylum"], useNA="ifany"), decreasing=TRUE))

cat("\n--- Reads per sample ---\n")
print(sort(rowSums(seqtab)))

cat("\nTotal ASVs:", ncol(seqtab), "\n")
cat("Samples:", nrow(seqtab), "\n")
