#!/usr/bin/env python
# coding: utf-8

import sys
import pandas as pd
import logging

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s: %(message)s')

GWAS, SNP_CONVERSIONS = sys.argv[1:]

# GWAS = '/net/topmed10/working/porchard/rnaseq/work/panukbb-gwas/results/biomarkers-30630-both_sexes-irnt.tsv.bgz'
# SNP_CONVERSIONS = '/net/topmed10/working/porchard/rnaseq/work/lift-panukbb-variants/results/conversion-table/snp-names.txt.gz'


logging.info(f'Loading {SNP_CONVERSIONS}')
snps = pd.read_csv(SNP_CONVERSIONS, sep='\t', compression='gzip')

logging.info(f'Loading {GWAS}')
gwas = pd.read_csv(GWAS, sep='\t', compression='gzip')
gwas['hg19_id'] = 'chr' + gwas.chr.astype(str) + '_' + gwas.pos.astype(str) + '_' + gwas.ref + '_' + gwas.alt
gwas = gwas.drop(columns=['chr', 'pos', 'ref', 'alt'])
gwas['missing'] = ~gwas.hg19_id.isin(set(snps.hg19_id))

PRESENT = (~gwas.missing).sum()
MISSING = gwas.missing.sum()
TOTAL = len(gwas)
assert(TOTAL == (PRESENT + MISSING))
logging.info('{:,} of {:,} SNPs are convertable'.format(PRESENT, TOTAL))

gwas = gwas.merge(snps[['hg19_id', 'hg38_id', 'chrom', 'start', 'end', 'ref', 'alt', 'swapped_alleles']])
assert(len(gwas) == PRESENT)

SWAP_SIGN_COLS = [i for i in gwas.columns if 'beta_' in i]
SUBTRACT_FROM_1_COLS = [i for i in gwas.columns if 'af_' in i]

logging.info('If ref and alt were swapped, will switch sign for these fields: {}'.format(', '.join(SWAP_SIGN_COLS)))
logging.info('If ref and alt were swapped, will recode x = 1 - x for these fields: {}'.format(', '.join(SUBTRACT_FROM_1_COLS)))

gwas.loc[gwas.swapped_alleles,SWAP_SIGN_COLS] = -1*gwas.loc[gwas.swapped_alleles,SWAP_SIGN_COLS]
gwas.loc[gwas.swapped_alleles,SUBTRACT_FROM_1_COLS] = 1-gwas.loc[gwas.swapped_alleles,SUBTRACT_FROM_1_COLS]


gwas = gwas.drop(columns=['hg19_id', 'hg38_id', 'swapped_alleles', 'missing'])
COL_ORDER = ['chrom', 'start', 'end', 'ref', 'alt']
OTHER_COLS = [i for i in gwas.columns if i not in COL_ORDER]
gwas = gwas[COL_ORDER + OTHER_COLS].rename(columns={'chrom': '#chrom'})

logging.info('Sorting')
gwas = gwas.sort_values(['#chrom', 'start'])

logging.info('Writing output')
gwas.astype(str).to_csv(sys.stdout, na_rep='NA', sep='\t', index=False)

logging.info('Done')