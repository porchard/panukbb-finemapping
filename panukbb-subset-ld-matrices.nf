#!/usr/bin/env nextflow

nextflow.enable.dsl=2

REGIONS = params.regions // chrom, start, end, ANCESTRY
LD_MAT_GLOB = params.ld_mat_glob // UKBB.{ANCESTRY}.ldadj.bm
LD_MAT_IDX_GLOB = params.ld_mat_idx_glob // UKBB.{ANCESTRY}.ldadj.variant.ht


process make_ld_mat {

    publishDir "${params.results}/ld-matrix"
    memory '10 GB'
    clusterOptions='--partition=topmed-working --exclude=topmed,topmed2,topmed[5-6],topmed[8-10]'
    maxForks 50
    maxRetries 2
    errorStrategy {task.attempt <= maxRetries ? 'retry' : 'ignore'}
    beforeScript 'source /net/topmed10/working/porchard/rnaseq/bin/config-pyspark.sh'
    cache 'lenient'

    input:
    tuple val(ancestry), path(ld_mat), path(ld_mat_idx), val(chrom), val(start), val(end)

    output:
    path("${chrom}_${start}_${end}___${ancestry}.ld.tsv.bgz")
    path("${chrom}_${start}_${end}___${ancestry}.variants.tsv.gz")

    """
    printf "${chrom}\\t${start}\\t${end}\\n" > regions.bed
    panukbb-subset-ld-matrix.py $ld_mat $ld_mat_idx regions.bed ''
    rm -rf *.bm
    mv ${chrom}_${start}_${end}.ld.tsv.bgz ${chrom}_${start}_${end}___${ancestry}.ld.tsv.bgz
    mv ${chrom}_${start}_${end}.variants.tsv.gz ${chrom}_${start}_${end}___${ancestry}.variants.tsv.gz
    """

}


workflow {
    ld_mat = Channel.fromPath(LD_MAT_GLOB, type: 'dir').map({it -> [it.getName().tokenize('.')[1], it]}) // ancestry, file
    ld_mat_idx = Channel.fromPath(LD_MAT_IDX_GLOB, type: 'dir').map({it -> [it.getName().tokenize('.')[1], it]}) // ancestry, file

    regions = Channel.from(file(REGIONS).readLines().collect({it -> it.replaceAll('\n', '').tokenize('\t')})).map({it -> [it[3], it[0], it[1], it[2]]}) // ancestry, chrom, start, end

    ld_mat.join(ld_mat_idx).combine(regions, by: 0) | make_ld_mat

}
