#!/usr/bin/env nextflow

nextflow.enable.dsl=2
URL_LIST = params.url_list

process dl {

    publishDir "${params.results}"
    maxForks 40
    memory '10 GB'

    input:
    val(url)

    output:
    path(outname)

    script:
    url_escape_pipe_character = url.replaceAll("\\|", "\\\\|")
    outname = url.replaceAll(".*/", "").replaceAll("\\|", '_')

    """
    aws s3 cp --only-show-errors --no-sign-request $url_escape_pipe_character $outname
    """

}


workflow {
    Channel.from(file(URL_LIST).readLines().collect({it -> it.replaceAll('\n', '')})) | dl
}