#!/usr/bin/env nextflow

nextflow.enable.dsl=2
GWAS_GLOB = params.gwas_glob
TRAITS = params.traits // TSV file listing traits to test for each ancestry. first column = trait_id, second column = ancestry (one line per trait and ancestry)
CHROM_SIZES = params.chrom_sizes

process get_regions_to_finemap {

    maxForks 40
    clusterOptions='--partition=topmed-working --exclude=topmed,topmed[2-9]'
    memory { 20.GB * task.attempt }
    maxRetries 3
    errorStrategy {task.attempt <= maxRetries ? 'retry' : 'ignore'}
    tag "${trait} ${ancestry}"
    container 'library://porchard/default/general:20220107'

    input:
    tuple val(trait), path(gwas), val(ancestry)

    output:
    tuple val(trait), val(ancestry), path("${trait}___${ancestry}.bed")

    """
    panukbb-get-regions-to-finemap-per-ancestry.py --ancestry $ancestry --gwas $gwas > ${trait}___${ancestry}.bed
    """

}


process merge_regions {

    maxForks 10
    executor 'local'
    cache 'lenient'
    tag "${trait} ${ancestry}"
    container 'library://porchard/default/general:20220107'

    input:
    tuple val(trait), val(ancestry), path('regions.bed')

    output:
    tuple val(trait), val(ancestry), path("${trait}___${ancestry}.regions.txt")

    """
    cat regions.bed | sort -k1,1 -k2n,2 | bedtools merge > ${trait}___${ancestry}.regions.txt
    """

}


process clip_regions {

    maxForks 10
    executor 'local'
    publishDir "${params.results}/regions-to-finemap-no-size-filter"
    cache 'lenient'
    tag "${trait} ${ancestry}"
    container 'library://porchard/default/general:20220107'

    input:
    tuple val(trait), val(ancestry), path('regions.bed'), path(chrom_sizes)

    output:
    tuple val(trait), val(ancestry), path("${trait}___${ancestry}.regions.txt")

    """
    chromSizesToBed $chrom_sizes | grep -v -e "random"  -e "chrUn" | perl -pe 's/^chr//' > chroms.bed
    bedtools intersect -a regions.bed -b chroms.bed > ${trait}___${ancestry}.regions.txt
    rm chroms.bed
    """

}


process size_filter {

    maxForks 10
    executor 'local'
    publishDir "${params.results}/regions-to-finemap-size-filter"
    cache 'lenient'
    tag "${trait} ${ancestry}"
    container 'library://porchard/default/general:20220107'

    input:
    tuple val(trait), val(ancestry), path('regions.bed')

    output:
    path("${trait}___${ancestry}.regions.txt")

    """
    cat regions.bed | awk '(\$3-\$2)<10000000' > ${trait}___${ancestry}.regions.txt
    """

}


workflow {
    ancestries_and_traits = Channel.from(file(TRAITS).readLines().collect({it -> it.replaceAll('\n', '').tokenize('\t')})) // trait_id, ancestry
    gwas = Channel.fromPath(GWAS_GLOB).map({it -> [it.getName().replaceAll('.tsv.bgz', ''), it]}) // trait_id, gwas
    merged = gwas.combine(ancestries_and_traits, by: 0) | get_regions_to_finemap | merge_regions
    merged.combine(Channel.fromPath(CHROM_SIZES)) | clip_regions | size_filter
}