#!/usr/bin/env nextflow

params.reads     = "$HOME/ngs_course/data/ecoli_reads/*_{1,2}.fastq.gz"
params.reference = "$HOME/ngs_course/data/practice/rel606.fna"
params.ploidy    = 1
params.outdir    = "results_variants"

process FASTP {
    tag "$id"
    publishDir "${params.outdir}/trimmed", mode: 'copy', pattern: "*.{html,json}"
    cpus 4

    input:
    tuple val(id), path(reads)

    output:
    tuple val(id), path("${id}_1.trim.fq.gz"), path("${id}_2.trim.fq.gz"), emit: trimmed
    path "${id}_fastp.json", emit: report

    script:
    """
    fastp -i ${reads[0]} -I ${reads[1]} \
          -o ${id}_1.trim.fq.gz -O ${id}_2.trim.fq.gz \
          --detect_adapter_for_pe --cut_tail --cut_tail_mean_quality 20 \
          --length_required 50 --thread ${task.cpus} \
          --html ${id}_fastp.html --json ${id}_fastp.json
    """
}

process BWA_INDEX {
    tag "index"

    input:
    path reference

    output:
    path "${reference}*"

    script:
    """
    bwa index ${reference}
    """
}

process BWA_MEM {
    tag "$id"
    cpus 4

    input:
    tuple val(id), path(r1), path(r2)
    path index
    val ref_name

    output:
    tuple val(id), path("${id}.sorted.bam")

    script:
    """
    bwa mem -t ${task.cpus} ${ref_name} ${r1} ${r2} \
      | samtools sort -@ 2 -o ${id}.sorted.bam -
    """
}

process MARKDUP {
    tag "$id"
    publishDir "${params.outdir}/aligned", mode: 'copy'
    cpus 4

    input:
    tuple val(id), path(bam)

    output:
    tuple val(id), path("${id}.markdup.bam"), path("${id}.markdup.bam.bai"), emit: bam
    path "${id}.flagstat.txt", emit: stats

    script:
    """
    samtools sort -n -@ 2 -o name.bam ${bam}
    samtools fixmate -m -@ 2 name.bam fix.bam
    samtools sort -@ 2 -o pos.bam fix.bam
    samtools markdup -@ 2 pos.bam ${id}.markdup.bam
    samtools index ${id}.markdup.bam
    samtools flagstat ${id}.markdup.bam > ${id}.flagstat.txt
    rm name.bam fix.bam pos.bam
    """
}

process CALL_VARIANTS {
    tag "$id"
    publishDir "${params.outdir}/variants", mode: 'copy'

    input:
    tuple val(id), path(bam), path(bai)
    path reference

    output:
    tuple val(id), path("${id}.filtered.vcf")

    script:
    """
    bcftools mpileup -Ou -f ${reference} --max-depth 250 \
                     --min-MQ 20 --min-BQ 20 ${bam} \
      | bcftools call -mv --ploidy ${params.ploidy} -Ov -o raw.vcf

    bcftools filter -i 'QUAL>=30 && DP>=10 && DP<=200' raw.vcf -o ${id}.filtered.vcf
    """
}

workflow {
    reads_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)
    ref_file = file(params.reference)
    ref_ch   = Channel.value(ref_file)

    trimmed = FASTP(reads_ch)
    index   = BWA_INDEX(ref_ch)

    aligned = BWA_MEM(trimmed.trimmed, index.collect(), ref_file.name)
    marked  = MARKDUP(aligned)
    vcfs    = CALL_VARIANTS(marked.bam, ref_ch)

    vcfs.view { id, vcf -> "Variants written: $id" }
}
