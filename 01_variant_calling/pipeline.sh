#!/usr/bin/env bash
# Bacterial variant calling: E. coli LTEE clone vs REL606 ancestor
# Usage: bash pipeline.sh
set -euo pipefail

REF=ref/rel606.fna
R1=data/SRR2584863_1.fastq.gz
R2=data/SRR2584863_2.fastq.gz
THREADS=4

mkdir -p results/{qc,trimmed,aligned,variants}

# 1. QC
fastqc "$R1" "$R2" -o results/qc -t "$THREADS"
multiqc results/qc -o results/qc

# 2. Trim
fastp -i "$R1" -I "$R2" \
      -o results/trimmed/R1.trim.fq.gz -O results/trimmed/R2.trim.fq.gz \
      --detect_adapter_for_pe --cut_tail --cut_tail_mean_quality 20 \
      --length_required 50 --thread "$THREADS" \
      --html results/trimmed/fastp.html --json results/trimmed/fastp.json

# 3. Align
bwa index "$REF"
bwa mem -t "$THREADS" "$REF" \
    results/trimmed/R1.trim.fq.gz results/trimmed/R2.trim.fq.gz \
  | samtools sort -@ 2 -o results/aligned/sorted.bam -
samtools index results/aligned/sorted.bam
samtools flagstat results/aligned/sorted.bam > results/aligned/flagstat.txt

# 4. Mark duplicates
samtools sort -n -@ 2 -o results/aligned/name.bam results/aligned/sorted.bam
samtools fixmate -m -@ 2 results/aligned/name.bam results/aligned/fix.bam
samtools sort -@ 2 -o results/aligned/pos.bam results/aligned/fix.bam
samtools markdup -@ 2 results/aligned/pos.bam results/aligned/markdup.bam
samtools index results/aligned/markdup.bam
rm results/aligned/{name,fix,pos}.bam

# 5. Call variants  -- NOTE: --ploidy 1, E. coli is haploid
bcftools mpileup -Ou -f "$REF" --max-depth 250 --min-MQ 20 --min-BQ 20 \
    results/aligned/markdup.bam \
  | bcftools call -mv --ploidy 1 -Ov -o results/variants/raw.vcf

# 6. Filter
bcftools filter -i 'QUAL>=30 && DP>=10 && DP<=200' \
    results/variants/raw.vcf -o results/variants/filtered.vcf

# 7. Annotate  -- NOTE: SnpEff DB names the sequence "Chromosome", VCF uses NC_012967.1
sed 's/^NC_012967.1/Chromosome/; s/ID=NC_012967.1/ID=Chromosome/' \
    results/variants/filtered.vcf > results/variants/renamed.vcf
snpEff -v Escherichia_coli_b_str_rel606 results/variants/renamed.vcf \
  > results/variants/annotated.vcf

echo "Done. Variants: $(grep -vc '^#' results/variants/annotated.vcf)"
