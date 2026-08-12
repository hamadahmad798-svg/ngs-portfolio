library(vegan)
seqtab <- readRDS("results/seqtab_nochim.rds")
taxa   <- readRDS("results/taxa.rds")

# Drop the mock - it's a control, not a biological sample
seqtab <- seqtab[rownames(seqtab) != "Mock", ]

# Mouse gut time course: Day 0-9 = early (pre-weaning), Day 141-150 = late
day <- as.numeric(sub("F3D", "", rownames(seqtab)))
group <- ifelse(day <= 9, "Early", "Late")
meta <- data.frame(sample=rownames(seqtab), day=day, group=group)
print(meta)

depth <- rowSums(seqtab)

cat("\n--- ALPHA DIVERSITY (raw, unrarefied) ---\n")
alpha <- data.frame(
  depth    = depth,
  observed = rowSums(seqtab > 0),
  shannon  = diversity(seqtab, "shannon"),
  group    = group
)
print(round(alpha[order(alpha$depth), c("depth","observed","shannon")], 3))

cat("\nCorrelation of observed richness with depth: ",
    round(cor(alpha$depth, alpha$observed, method="spearman"), 3), "\n")
cat("Correlation of Shannon with depth:            ",
    round(cor(alpha$depth, alpha$shannon,  method="spearman"), 3), "\n")

cat("\n--- ALPHA DIVERSITY (rarefied to min depth) ---\n")
set.seed(42)
rare <- rrarefy(seqtab, min(depth))
alpha$obs_rare <- rowSums(rare > 0)
alpha$sha_rare <- diversity(rare, "shannon")
print(round(alpha[, c("depth","observed","obs_rare","shannon","sha_rare")], 3))

cat("\n--- Early vs Late (rarefied) ---\n")
print(t.test(obs_rare ~ group, data=alpha))
print(t.test(sha_rare ~ group, data=alpha))

cat("\n--- BETA DIVERSITY ---\n")
rel <- sweep(seqtab, 1, rowSums(seqtab), "/")
bc  <- vegdist(rel, method="bray")

set.seed(42)
print(adonis2(bc ~ group, data=meta, permutations=999))

nmds <- metaMDS(bc, k=2, trace=FALSE)
cat("\nNMDS stress:", round(nmds$stress, 4), "\n")

png("results/nmds.png", width=800, height=600)
cols <- ifelse(group=="Early", "#D85A30", "#2C7A7B")
plot(nmds$points, col=cols, pch=19, cex=2,
     xlab="NMDS1", ylab="NMDS2", main="Bray-Curtis NMDS")
text(nmds$points, labels=rownames(seqtab), pos=3, cex=0.7)
legend("topright", legend=c("Early","Late"), col=c("#D85A30","#2C7A7B"), pch=19)
dev.off()

saveRDS(list(meta=meta, alpha=alpha, bc=bc), "results/diversity.rds")
EOFcat > diversity.R << 'EOF'
library(vegan)
seqtab <- readRDS("results/seqtab_nochim.rds")
taxa   <- readRDS("results/taxa.rds")

# Drop the mock - it's a control, not a biological sample
seqtab <- seqtab[rownames(seqtab) != "Mock", ]

# Mouse gut time course: Day 0-9 = early (pre-weaning), Day 141-150 = late
day <- as.numeric(sub("F3D", "", rownames(seqtab)))
group <- ifelse(day <= 9, "Early", "Late")
meta <- data.frame(sample=rownames(seqtab), day=day, group=group)
print(meta)

depth <- rowSums(seqtab)

cat("\n--- ALPHA DIVERSITY (raw, unrarefied) ---\n")
alpha <- data.frame(
  depth    = depth,
  observed = rowSums(seqtab > 0),
  shannon  = diversity(seqtab, "shannon"),
  group    = group
)
print(round(alpha[order(alpha$depth), c("depth","observed","shannon")], 3))

cat("\nCorrelation of observed richness with depth: ",
    round(cor(alpha$depth, alpha$observed, method="spearman"), 3), "\n")
cat("Correlation of Shannon with depth:            ",
    round(cor(alpha$depth, alpha$shannon,  method="spearman"), 3), "\n")

cat("\n--- ALPHA DIVERSITY (rarefied to min depth) ---\n")
set.seed(42)
rare <- rrarefy(seqtab, min(depth))
alpha$obs_rare <- rowSums(rare > 0)
alpha$sha_rare <- diversity(rare, "shannon")
print(round(alpha[, c("depth","observed","obs_rare","shannon","sha_rare")], 3))

cat("\n--- Early vs Late (rarefied) ---\n")
print(t.test(obs_rare ~ group, data=alpha))
print(t.test(sha_rare ~ group, data=alpha))

cat("\n--- BETA DIVERSITY ---\n")
rel <- sweep(seqtab, 1, rowSums(seqtab), "/")
bc  <- vegdist(rel, method="bray")

set.seed(42)
print(adonis2(bc ~ group, data=meta, permutations=999))

nmds <- metaMDS(bc, k=2, trace=FALSE)
cat("\nNMDS stress:", round(nmds$stress, 4), "\n")

png("results/nmds.png", width=800, height=600)
cols <- ifelse(group=="Early", "#D85A30", "#2C7A7B")
plot(nmds$points, col=cols, pch=19, cex=2,
     xlab="NMDS1", ylab="NMDS2", main="Bray-Curtis NMDS")
text(nmds$points, labels=rownames(seqtab), pos=3, cex=0.7)
legend("topright", legend=c("Early","Late"), col=c("#D85A30","#2C7A7B"), pch=19)
dev.off()

saveRDS(list(meta=meta, alpha=alpha, bc=bc), "results/diversity.rds")
