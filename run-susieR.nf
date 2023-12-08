#!/usr/bin/env nextflow

nextflow.enable.dsl=2

TO_RUN = params.to_run
SNPS_EXCLUDE = params.snps_exclude

process run_susie {

    publishDir "${params.results}/susie-out"
    //memory { 20.GB * task.attempt }
    //memory { 20.GB + (100.GB * task.attempt - 1) }
    memory { 120.GB * task.attempt }
    //clusterOptions='--partition=topmed-working --exclude=topmed,topmed[2-10]'
    clusterOptions='--partition=main'
    maxForks 50
    container 'docker.io/porchard/coloc:20220524'
    maxRetries 4
    time '14d'
    errorStrategy {task.attempt <= maxRetries ? 'retry' : 'ignore'}

    input:
    tuple val(trait), path(gwas), path(gwas_idx), val(ancestry), val(chrom), val(start), val(end), path(ld_matrix), path(ld_variants), path(snps_exclude)

    output:
    path("${prefix}.rda")

    script:
    prefix = trait + '___' + ancestry + '___' +  chrom + '_' + start + '_' + end

    """
    run-susieR.R $gwas $ancestry $ld_matrix $ld_variants ${chrom}:${start}-${end} $snps_exclude ${prefix}.rda ${prefix}.png
    """

}



workflow {
    to_run = Channel.fromPath(TO_RUN)
    susie_out = Channel.from(file(TO_RUN).readLines().collect({it -> it.replaceAll('\n', '').tokenize('\t')})).combine(Channel.fromPath(SNPS_EXCLUDE)) | run_susie
}