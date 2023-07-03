#!/usr/bin/env nextflow

nextflow.enable.dsl=2
ANCESTRIES = ['EUR', 'AFR']

process fetch_matrix {

    publishDir "${params.results}/full-matrix"
    memory '10 GB'
    errorStrategy 'ignore'
    time '168h'

    input:
    val(ancestry)

    output:
    tuple val(ancestry), path("UKBB.${ancestry}.ldadj.bm")
    
    """
    aws s3 sync --only-show-errors --no-sign-request s3://pan-ukb-us-east-1/ld_release/UKBB.${ancestry}.ldadj.bm UKBB.${ancestry}.ldadj.bm
    """

}


process fetch_index {

    publishDir "${params.results}/full-index"
    memory '10 GB'
    errorStrategy 'ignore'

    input:
    val(ancestry)

    output:
    tuple val(ancestry), path("UKBB.${ancestry}.ldadj.variant.ht")
    
    """
    aws s3 sync --only-show-errors --no-sign-request s3://pan-ukb-us-east-1/ld_release/UKBB.${ancestry}.ldadj.variant.ht UKBB.${ancestry}.ldadj.variant.ht
    """

}


workflow {
    a = Channel.from(ANCESTRIES)
    fetch_matrix(a)
    fetch_index(a)
}