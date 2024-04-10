#!/usr/bin/env nextflow

nextflow.enable.dsl=2

SUSIE_GLOB = params.susie_glob


process credible_sets_and_convergence {

    container 'docker.io/porchard/coloc:20220505'
    clusterOptions '--partition=topmed-working --exclude=topmed,topmed[2-10]'
    maxForks 5

    input:
    path("rda/*")

    output:
    path("gwas.cs.txt"), emit: cs
    path("gwas.converged.txt"), emit: converged

    """
    ls rda/* > rda-files.txt
    susieR-get-credible-sets-and-convergence.R rda-files.txt gwas.
    """

}


process concat_cs_and_convergence {

    publishDir "${params.results}/susie-cs-and-convergence"
    clusterOptions '--partition=topmed-working --exclude=topmed,topmed[2-10]'
    container 'library://porchard/default/general:20220107'

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



workflow {
    susie_out = Channel.fromPath(SUSIE_GLOB)

    before_concat = susie_out.collate(500) | credible_sets_and_convergence

    concat_cs_and_convergence(before_concat.cs.toSortedList(), before_concat.converged.toSortedList())

}
