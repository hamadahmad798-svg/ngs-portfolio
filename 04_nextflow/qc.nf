#!/usr/bin/env nextflow

params.reads   = "$HOME/ngs_course/microbiome/data/MiSeq_SOP/*_R{1,2}_001.fastq"
params.outdir  = "results"

process FASTQC {
    tag "$sample_id"
    publishDir "${params.outdir}/fastqc", mode: 'copy'
    cpus 2

    input:
    tuple val(sample_id), path(reads)

    output:
    path "*_fastqc.{zip,html}"

    script:
    """
    fastqc ${reads} -t ${task.cpus}
    """
}

process MULTIQC {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path '*'

    output:
    path "multiqc_report.html"

    script:
    """
    multiqc .
    """
}

workflow {
    read_pairs = Channel.fromFilePairs(params.reads, checkIfExists: true)
    read_pairs.view { id, files -> "Found: $id" }

    fastqc_out = FASTQC(read_pairs)
    MULTIQC(fastqc_out.collect())
}
