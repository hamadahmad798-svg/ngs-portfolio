seqtab <- readRDS("results/seqtab_nochim.rds")
taxa   <- readRDS("results/taxa.rds")
seqtab <- seqtab[rownames(seqtab) != "Mock", ]

day <- as.numeric(sub("F3D", "", rownames(seqtab)))
group <- ifelse(day <= 9, "Early", "Late")

# Collapse ASVs to genus, then convert to relative abundance
gen <- taxa[colnames(seqtab), "Genus"]
gen[is.na(gen)] <- "Unassigned"
agg <- t(rowsum(t(seqtab), group = gen))
rel <- sweep(agg, 1, rowSums(agg), "/") * 100

mE <- colMeans(rel[group == "Early", ])
mL <- colMeans(rel[group == "Late",  ])

df <- data.frame(genus = names(mE), Early = round(mE,2), Late = round(mL,2),
                 diff = round(mL - mE, 2))
df <- df[order(-abs(df$diff)), ]

cat("--- Top 20 genera by change in mean relative abundance (%) ---\n")
print(head(df, 20), row.names = FALSE)

cat("\n--- Wilcoxon tests, genera with mean abundance > 0.5% ---\n")
keep <- names(which(colMeans(rel) > 0.5))
res <- data.frame()
for (g in keep) {
  p <- suppressWarnings(wilcox.test(rel[,g] ~ group)$p.value)
  res <- rbind(res, data.frame(genus=g, Early=round(mE[g],2),
                               Late=round(mL[g],2), p=signif(p,3)))
}
res$padj <- signif(p.adjust(res$p, method="BH"), 3)
print(res[order(res$p), ], row.names = FALSE)
