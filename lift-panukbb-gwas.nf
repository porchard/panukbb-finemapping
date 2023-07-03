#!/usr/bin/env nextflow

nextflow.enable.dsl=2

GWAS_GLOB = params.gwas_glob
SNP_CONVERSIONS = params.snp_conversions

process lift {

    container 'library://porchard/default/general:20220107'
    publishDir "${params.results}/lifted"
    clusterOptions='--partition=topmed-working --exclude=topmed,topmed2,topmed[5-6],topmed[8-10]'
    maxForks 10
    memory { 50.GB * task.attempt }
    maxRetries 2
    errorStrategy {task.attempt <= maxRetries ? 'retry' : 'ignore'}

    input:
    each path(gwas)
    path(snp_conversions)

    output:
    path(outfile)

    script:
    outfile = gwas.getName().replaceAll('.tsv.bgz', '.hg38.tsv.bgz')

    """
    lift-panukbb-gwas.py $gwas $snp_conversions | bgzip -c > $outfile
    """

}


process index {

    container 'library://porchard/default/general:20220107'
    publishDir "${params.results}/lifted"
    clusterOptions='--partition=topmed-working --exclude=topmed,topmed2,topmed[5-6],topmed[8-10]'
    maxForks 30

    input:
    path(x)

    output:
    path("*.tbi")

    """
    tabix -p bed $x
    """

}

workflow {
    lift(Channel.fromPath(GWAS_GLOB), Channel.fromPath(SNP_CONVERSIONS)) | index
}