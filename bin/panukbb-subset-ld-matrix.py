#!/usr/bin/env python
# coding: utf-8

import os
import sys
import subprocess
import hail as hl
from hail.linalg import BlockMatrix
import pandas as pd
import gzip

LD_MATRIX, VARIANT_INDEX, REGIONS, PREFIX = sys.argv[1:]

#LD_MATRIX = '/net/topmed10/working/porchard/rnaseq/work/panukbb-ld-matrices/results/full-matrix/UKBB.EUR.ldadj.bm'
#VARIANT_INDEX = '/net/topmed10/working/porchard/rnaseq/work/panukbb-ld-matrices/results/full-index/UKBB.EUR.ldadj.variant.ht'
#GWAS = '/net/topmed10/working/porchard/rnaseq/data/panukbb/gwas/continuous-30140-both_sexes-irnt.tsv.bgz'
#ANCESTRY = 'EUR'
#PREFIX = 'test-susie-in'

bm = BlockMatrix.read(LD_MATRIX)
ht_idx = hl.read_table(VARIANT_INDEX)

def subset_ld_matrix_to_locus(ldmat, ldmat_idx, region):
    """"
    Region arg is like bed, i.e. 0 indexed w/ start = inclusive and end = exclusive. Follows format {chr}:{start}-{end}
    """
    # https://hail.is/docs/0.2/functions/genetics.html#hail.expr.functions.parse_locus_interval
    # parse_locus_interval has min start = 1 when left is inclusive, so it's indexed from 1 (unlike bed which indexes from 0)
    # therefore, add 1 to start and end
    chrom = region.split(':')[0]
    start, end = region.split(':')[1].split('-')
    start = int(start) + 1
    end = int(end) + 1
    new_region = f'{chrom}:{start}-{end}'

    interval = hl.parse_locus_interval(new_region)
    ldmat_idx_subset = ldmat_idx.filter(interval.contains(ldmat_idx.locus))
    variants_df = ldmat_idx_subset.to_pandas()
    variants_df['variant_id'] = variants_df.locus.astype(str).str.replace(':', '_') + '_' + variants_df.alleles.map(lambda x: '_'.join(x))
    idx = variants_df.idx.to_list()
    ldmat_subset = ldmat.filter(idx, idx)
    return (ldmat_subset, variants_df)


regions = pd.read_csv(REGIONS, header=None, sep='\t', names=['chrom', 'start', 'end'])
regions['region'] = regions.chrom.astype(str) + ':' + regions.start.astype(str) + '-' + regions.end.astype(str)

for region in regions.region.to_list():
    (ld, df) = subset_ld_matrix_to_locus(bm, ht_idx, region)
    # export the LD matrix and GWAS subset
    region_name = region.replace(':', '_').replace('-', '_')
    BM_EXPORT = f'{PREFIX}{region_name}.bm'
    TXT_EXPORT = f'{PREFIX}{region_name}.ld.tsv.bgz'
    DF_EXPORT = f'{PREFIX}{region_name}.variants.tsv.gz'

    ld.write(BM_EXPORT, force_row_major=True)
    BlockMatrix.export(
        BM_EXPORT,
        TXT_EXPORT,
        delimiter='\t'
    )

    df.to_csv(DF_EXPORT, sep='\t', index=False, compression='gzip')
