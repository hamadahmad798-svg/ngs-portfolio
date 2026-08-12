#!/usr/bin/env nextflow

params.greeting = "Hello"

process SAY {
    input:
    val name

    output:
    stdout

    script:
    """
    echo "${params.greeting}, ${name}!"
    """
}

workflow {
    names = Channel.of('Hamad', 'REL606', 'MazF')
    SAY(names) | view
}
