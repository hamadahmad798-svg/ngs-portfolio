library(dada2)
seqtab <- readRDS("results/seqtab_nochim.rds")

taxa <- assignTaxonomy(seqtab, "ref/silva_nr99_v138.1_train_set.fa.gz", multithread=4)
saveRDS(taxa, "results/taxa.rds")

cat("\n--- Taxonomic assignment success by rank ---\n")
print(colSums(!is.na(taxa)))
cat("Total ASVs:", nrow(taxa), "\n")

cat("\n--- Phylum distribution ---\n")
print(sort(table(taxa[,"Phylum"], useNA="ifany"), decreasing=TRUE))

# Mock community check
mock <- seqtab["Mock",]
mock <- mock[mock > 0]
cat("\n--- Mock community ---\n")
cat("ASVs detected in Mock:", length(mock), "\n")
cat("Reads in Mock:", sum(mock), "\n")
cat("Abundance of each ASV in Mock:\n")
print(sort(mock, decreasing=TRUE))

mock.tax <- taxa[names(mock), c("Family","Genus")]
print(mock.tax)
