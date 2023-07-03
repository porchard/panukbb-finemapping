#!/usr/bin/env nextflow

nextflow.enable.dsl=2

SUSIE_GLOB = params.susie_glob
TO_RUN = params.to_run
CHAIN = params.chain
SNP_CONVERSIONS = params.snp_conversions
SNP_CONVERSIONS_IDX = params.snp_conversions_idx


process lift_regions {

    container 'library://porchard/default/general:20220107'
    publishDir "${params.results}/lifted-regions"
    executor 'local'

    input:
    path(to_run)
    path(chain)

    output:
    path('regions-hg38.bed')

    """
    cut -f5-7 $to_run | perl -pe 's/^/chr/' | sort -k1,1 -k2n,2 > tmp.txt
    cat tmp.txt | perl -pe 's/\\t/:/g' | paste tmp.txt - > regions-hg19.bed
    liftOver regions-hg19.bed $chain regions-hg38.bed unmapped.bed
    """

}


process lift_susie {

    publishDir "${params.results}/susie-hg38"
    memory '10 GB'
    //queue 'nomosix'
    container 'docker.io/porchard/coloc:20220505'
    //executor 'local'
    clusterOptions='--partition=topmed-working --exclude=topmed,topmed[2-10]'

    input:
    each path("rda/*")
    path(lifted_regions)
    path(snp_conversions)
    path(snp_conversions_idx)

    output:
    path("*.rda"), emit: for_cs
    path("susie-conversion-summary.txt"), emit: summary

    """
    ls rda/* > rda-files.txt
    lift-susie.R $snp_conversions $lifted_regions rda-files.txt
    """

}


process credible_sets_and_convergence {

    container 'docker.io/porchard/coloc:20220505'
    clusterOptions='--partition=topmed-working --exclude=topmed,topmed[2-10]'

    input:
    path("rda/*")

    output:
    path("gwas.cs.txt"), emit: cs
    path("gwas.converged.txt"), emit: converged

    """
    ls rda/* > rda-files.txt
    susieR-get-credible-sets-and-convergence-2.R rda-files.txt gwas.
    """

}


process concat_summary {

    publishDir "${params.results}/susie-liftover-summary"
    executor 'local'

    input:
    path("summary.*.txt")

    output:
    path('summary.txt')


    """
    cat summary.1.txt | awk 'NR==1' > header.txt
    cp header.txt summary.txt
    cat summary.*.txt | grep -v -f header.txt >> summary.txt 
    """

}

process concat_cs_and_convergence {

    publishDir "${params.results}/susie-cs-and-convergence"
    executor 'local'

    input:
    path("cs.*.txt")
    path("converged.*.txt")

    output:
    path('gwas.cs.txt')
    path('gwas.converged.txt')


    """
    cat cs.*.txt | perl -pe 's/.rda//' > gwas.cs.txt
    cat converged.*.txt | perl -pe 's/.rda//' > gwas.converged.txt 
    """

}


process lbf {

    executor 'local'
    container 'docker.io/porchard/coloc:20220505'
    maxForks 10
    memory '30 GB'

    input:
    tuple val(trait), path(rda)

    output:
    path("${trait}.bed")

    """
    susieR-to-lbf-bed.R $trait ${trait}.bed ${rda.join(' ')}
    """

}


process concat_lbf {

    publishDir "${params.results}/susie-lbf"
    container 'library://porchard/default/general:20220107'
    clusterOptions='--partition=topmed-working --exclude=topmed,topmed2,topmed[4-10]'
    memory '100 GB'

    input:
    path(x)

    output:
    path('lbf.bed.gz')
    path('lbf.bed.gz.tbi')

    """
    merge-lbf-tabix.R lbf.tmp ${x.join(' ')}
    cat lbf.tmp | awk 'NR==1' > lbf.bed
    cat lbf.tmp | grep -v '^#' | sort -k1,1 -k2n,2 >> lbf.bed
    bgzip lbf.bed
    tabix -p bed lbf.bed.gz
    """

}


workflow {
    to_run = Channel.fromPath(TO_RUN)
    chain = Channel.fromPath(CHAIN)
    snp_conversions = Channel.fromPath(SNP_CONVERSIONS)
    snp_conversions_idx = Channel.fromPath(SNP_CONVERSIONS_IDX)
    susie_out = Channel.fromPath(SUSIE_GLOB)
    
    lifted_regions = lift_regions(to_run, chain)
    lifted = lift_susie(susie_out.collate(500), lifted_regions, snp_conversions, snp_conversions_idx)
    lifted.summary.toSortedList() | concat_summary
    before_concat = credible_sets_and_convergence(lifted.for_cs.flatten().collate(500))
    concat_cs_and_convergence(before_concat.cs.toSortedList(), before_concat.converged.toSortedList())


    (lifted.for_cs.flatten().map({it -> [it.getName().split('___')[0], it]}).groupTuple() | lbf).toSortedList() | concat_lbf

}
