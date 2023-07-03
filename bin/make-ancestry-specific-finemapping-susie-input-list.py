#!/usr/bin/env python
# coding: utf-8

import os
import re
import sys
import glob
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--gwas-glob', dest='gwas_glob', help='Glob for GWAS files. Index files are assumed to be in the same location')
parser.add_argument('--regions-glob', dest='regions_glob', help='Glob for files indicating which regions to finemap for each GWAS. Filenames formatted as {trait}___{ancestry}.regions.txt')
parser.add_argument('--ld-matrix-glob', dest='ld_matrix_glob', nargs='+', help='')
parser.add_argument('--ld-matrix-variants-glob', dest='ld_matrix_variants_glob', nargs='+', help='')
args = parser.parse_args()


GWAS_GLOB = args.gwas_glob
REGIONS_GLOB = args.regions_glob
LD_MATRIX_GLOB = args.ld_matrix_glob
LD_MATRIX_VARIANTS_GLOB = args.ld_matrix_variants_glob


# GWAS_GLOB = '/net/topmed10/working/porchard/rnaseq/work/panukbb-gwas/results/*.bgz'
# ##GWAS_INDEX_GLOB = '/net/topmed10/working/porchard/rnaseq/work/panukbb-gwas/results/*.tbi'
# REGIONS_GLOB = '/net/topmed10/working/porchard/rnaseq/work/panukbb-extract-ld-matrices/results/merge-regions/*'
# LD_MATRIX_GLOB = '/net/topmed10/working/porchard/rnaseq/work/panukbb-extract-ld-matrices/results/ld-matrix/*.ld.tsv.bgz'
# LD_MATRIX_VARIANTS_GLOB = '/net/topmed10/working/porchard/rnaseq/work/panukbb-extract-ld-matrices/results/ld-matrix/*.variants.tsv.gz'


# make table:
# gwas
# gwas index
# chrom
# start
# end
# ld matrix
# ld variants

def parse_ld_filename(f):
    # '{chrom}_{start}_{end}___{ancestry}.{ld.tsv.bgz,variants.tsv.gz}
    RE = re.compile('^(.*)_(.*)_(.*)___(.*).(ld.tsv.bgz|variants.tsv.gz)$')
    chrom, start, end, ancestry, _ = RE.match(os.path.basename(f)).groups()
    return {'chrom': chrom, 'start': start, 'end': end, 'ancestry': ancestry}

gwas = pd.DataFrame({'gwas': glob.glob(GWAS_GLOB)})
gwas['gwas_idx'] = gwas.gwas + '.tbi'
gwas['trait'] = gwas.gwas.map(lambda f: os.path.basename(f).replace('.tsv.bgz', ''))

regions = pd.concat([pd.read_csv(f, sep='\t', header=None, names=['chrom', 'start', 'end'], dtype=str).assign(trait=os.path.basename(f).replace('.regions.txt', '').split('___')[0], ancestry=os.path.basename(f).replace('.regions.txt', '').split('___')[1]) for f in glob.glob(REGIONS_GLOB)])

assert(all(regions.trait.isin(gwas.trait.to_list())))

ld_matrices = []
for i in LD_MATRIX_GLOB:
    ld_matrices += glob.glob(i)
ld_variants = []
for i in LD_MATRIX_VARIANTS_GLOB:
    ld_variants += glob.glob(i)
ld_matrices = pd.DataFrame({'matrix': ld_matrices})
ld_variants = pd.DataFrame({'variants': ld_variants})
for i in ['chrom', 'start', 'end', 'ancestry']:
    ld_matrices[i] = ld_matrices.matrix.map(lambda x: parse_ld_filename(x)[i])
    ld_variants[i] = ld_variants.variants.map(lambda x: parse_ld_filename(x)[i])
assert(len(ld_matrices) == len(ld_variants))

ld = ld_matrices.merge(ld_variants)
run = gwas[['trait', 'gwas', 'gwas_idx']].merge(regions).merge(ld)
run[['trait', 'gwas', 'gwas_idx', 'ancestry', 'chrom', 'start', 'end', 'matrix', 'variants']].to_csv(sys.stdout, sep='\t', index=False)
