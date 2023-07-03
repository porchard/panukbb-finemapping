#!/usr/bin/env nextflow

nextflow.enable.dsl=2

VARIANT_FILE = params.variant_file
CHAIN_FILE = params.chain_file
FASTA = params.fasta

process make_vcf {

    publishDir "${params.results}/vcf"
    container 'library://porchard/default/general:20220107'

    input:
    path(x)

    output:
    path('panukbb-hg19.vcf')

    """
    #!/usr/bin/env python
    # coding: utf-8

    import csv
    import gzip

    fh_out = open('panukbb-hg19.vcf', 'w')

    fh_out.write("##fileformat=VCFv4.2\\n")
    fh_out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\\n')
    fh_out.write('\\t'.join(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'fakesample']) + '\\n')

    with gzip.open('${x}', 'rt') as f:
        reader = csv.DictReader(f, delimiter='\\t')
        for line in reader:
            fh_out.write('\\t'.join(['chr' + line['chrom'], line['pos'], line['varid'], line['ref'], line['alt'], '.', 'PASS', '.', 'GT', './.']) + '\\n')
    """

}

process make_sequence_dict {

    container 'library://porchard/default/general:20220107'
    memory '8 GB'

    input:
    path('hg38.fa')

    output:
    tuple path('hg38.fa'), path('hg38.dict')

    """
    java -Xmx4g -Xms4g -jar \$PICARD_JAR CreateSequenceDictionary R=hg38.fa O=hg38.dict
    """

}


process lift_vcf {

    publishDir "${params.results}/vcf"
    container 'library://porchard/default/general:20220107'
    memory '25 GB'
    time '168h'

    input:
    path(vcf)
    path(chain)
    tuple path(fasta), path(fasta_dict)

    output:
    path('panukbb-hg38.vcf.gz'), emit: lifted
    path('rejected.vcf.gz'), emit: rejected

    """
    java -Xmx20g -Xms20g -jar \$PICARD_JAR LiftoverVcf --CHAIN $chain --INPUT $vcf --OUTPUT panukbb-hg38.vcf.gz --REJECT rejected.vcf.gz --WARN_ON_MISSING_CONTIG true --REFERENCE_SEQUENCE $fasta --RECOVER_SWAPPED_REF_ALT true --WRITE_ORIGINAL_ALLELES true --WRITE_ORIGINAL_POSITION true
    """

}

process make_conversions {

    publishDir "${params.results}/conversion-table"
    memory '25 GB'
    container 'docker.io/porchard/general:20220406125608'

    input:
    path(vcf)

    output:
    path("snp-names.txt.gz")
    path("snp-names.txt.gz.tbi")

    """
    bcftools view --drop-genotypes --no-header $vcf | cut -f1-5,8 > variant-info.txt
    make-snp-conversions.py variant-info.txt > snp-names.txt
    bgzip snp-names.txt
    tabix -p bed snp-names.txt.gz
    """

}


workflow {
    vcf = Channel.fromPath(VARIANT_FILE) | make_vcf
    chain = Channel.fromPath(CHAIN_FILE)
    fasta = Channel.fromPath(FASTA)

    lift_vcf(vcf, chain, make_sequence_dict(fasta)).lifted | make_conversions
}