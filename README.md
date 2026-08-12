# NGS Bioinformatics Portfolio

Reproducible analyses of public sequencing data: bacterial variant calling, RNA-seq
differential expression, and 16S microbiome profiling. Each includes the diagnostic
checks used to separate real signal from methodological artefact.

**Hamad Ahmad** — BS Biotechnology, Capital University of Science and Technology, Islamabad

---

## 1. Bacterial variant calling — *E. coli* LTEE clone

**Data:** SRR2584863 · **Reference:** *E. coli* B str. REL606 (GCF_000017985.1)

**Results:** 26 high-confidence variants (20 SNPs, 6 indels), 83x mean depth,
99.993% of genome at >=10x. Three frameshifts, 17 missense, **zero synonymous** —
consistent with positive selection, though 20 sites is too few for a formal dN/dS test.

Affected genes are enriched for global regulators (*topA*, *fis*, *malT*, *iclR*,
*glpR*, *nadR*) rather than structural enzymes — adaptation via regulatory rewiring.

**Reference validation:**

| Metric | REL606 (correct) | K-12 (incorrect) |
|---|---|---|
| Mapped | 99.996% | 94.49% |
| Properly paired | 98.17% | 91.25% |
| Singletons | 73 | 19,326 |

Mapping rate fell 5.5 pp; singletons rose 265-fold. Reference suitability should be
judged on paired-end statistics, not mapping rate — 94.5% would pass most QC thresholds.

**Errors caught:** default diploid model produced 5 spurious calls in a haploid organism,
with no warning; a sequence-name mismatch silently annotated all variants as intergenic.

`01_variant_calling/`

---

## 2. RNA-seq differential expression — MazF toxin induction

**Data:** PRJDB5742, *E. coli* K-12 · 3 control vs 3 mazF-induced (20 min)

**Internal control:** *mazF* itself at log2FC +7.09, confirming correct labelling
and successful induction.

**Normalisation assumption tested and violated:** 1,494 of 2,663 transcripts reached
padj < 0.05. Median-of-ratios normalisation assumes most features are unchanged; at 56%
significant that fails, so fold-change magnitudes are not quantitative.

**Two effects separated:**

| Expression quartile | median log2FC | ACA density (/kb) |
|---|---|---|
| Q1 (lowest) | +1.36 | 10.98 |
| Q2 | +0.10 | 10.79 |
| Q3 | -0.11 | 10.54 |
| Q4 (highest) | -0.68 | 10.87 |

Fold-change declines monotonically with expression — compositional re-centring under
global mRNA loss. ACA density (MazF's cleavage motif) is flat across quartiles, so the
enrichment in down-regulated transcripts (11.83 vs 10.53/kb) is not confounded by
expression and supports genuine sequence-specific cleavage.

**Limitations:** underpowered (165k-427k mapped reads/sample); library size confounded
with condition; one control replicate is a PCA outlier not explained by depth.

`02_rnaseq/`

---

## 3. 16S microbiome — mouse gut weaning transition

**Data:** Schloss MiSeq SOP · 19 samples + mock community · DADA2 + SILVA

**Validation:** the mock community returned the expected strain count with zero chimeras.
Three *Streptococcus* and two *Bacteroides* species resolved as five distinct ASVs —
a distinction 97% OTU clustering would have collapsed.

**Results:** 232 ASVs; 3.6% of reads chimeric (20.8% of ASVs — numerous but rare).
Alpha diversity unchanged (Shannon p = 0.146) while composition differed strongly
(PERMANOVA R2 = 0.43, p = 0.001; NMDS stress 0.080). The community was replaced,
not depleted.

**Dispersion check:** betadisper p = 0.006 (0.224 Early vs 0.153 Late) — pre-weaning
communities are more heterogeneous. This qualifies the PERMANOVA and is itself a finding.

**Resolution:** 100% to Order, 90.5% to Family, 56% to Genus, 0% to Species. The
unassigned genus fraction (60% of reads) was almost entirely Muribaculaceae. Results
reported at family level.

`03_microbiome/`

---

## 4. Reproducible pipelines — Nextflow + Docker

- `qc.nf` — FastQC + MultiQC, arbitrary sample count
- `variants.nf` — fastp -> BWA -> markdup -> bcftools

**Validated against manual analysis:** identical variant calls (`diff` on chrom, pos,
ref, alt returned nothing). A 4-read difference in 2.8M from fastp version affected
no call.

`04_nextflow/`

---

## Tools

**Environment:** Linux · Bash · Conda/Bioconda · Git · Nextflow · Docker

**NGS:** FastQC · MultiQC · fastp · BWA-MEM · SAMtools · BCFtools · SnpEff · IGV · Salmon · seqkit

**R:** DESeq2 · tximport · DADA2 · vegan

**Python:** pandas · SciPy

## Approach

Every analysis here includes a component where the answer is known in advance — a mock
community, an induced gene, a validated reference — and every result states the
conditions under which it would not hold. Sequencing measures relative quantities;
distinguishing biological change from compositional artefact is usually the analysis,
not a footnote.

## Notes

Raw data not included (see `.gitignore`); all datasets are public with accessions above.
Analyses are exploratory and reproducible, not clinical.
