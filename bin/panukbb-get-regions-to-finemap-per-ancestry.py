#!/usr/bin/env python
# coding: utf-8

import sys
import pandas as pd
import numpy as np
import logging
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--extension', default=250000, type=int)
parser.add_argument('--ancestry', required=True, choices=['AFR', 'AMR', 'CSA', 'EAS', 'EUR', 'MID'])
parser.add_argument('--gwas', required=True)
args = parser.parse_args()

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s: %(message)s')


EXTENSION = args.extension
GWAS = args.gwas
ancestry = args.ancestry

logging.info('Loading gwas')
gwas = pd.read_csv(GWAS, sep='\t', compression='gzip', usecols=['chr', 'pos', 'ref', 'alt', f'beta_{ancestry}', f'pval_{ancestry}', f'low_confidence_{ancestry}'], dtype={'chr': 'str', 'pos': 'int', 'ref': 'str', 'alt': 'str', f'beta_{ancestry}': 'float', f'pval_{ancestry}': 'float', f'low_confidence_{ancestry}': 'str'})
gwas_ancestry = gwas.rename(columns=lambda x: x.replace(f'_{ancestry}', ''))
logging.info('Finished loading gwas')

if 'beta' not in gwas_ancestry or (~gwas_ancestry.beta.isnull()).sum() == 0:
    raise ValueError(f'No GWAS stats for ancestry {ancestry}')

gwas_ancestry = gwas_ancestry[~gwas_ancestry.beta.isnull()]
gwas_ancestry = gwas_ancestry[gwas_ancestry.low_confidence=='false']

# get genome-wide significant
sig = gwas_ancestry[gwas_ancestry.pval<=np.log(5e-8)].sort_values('pval')
sig['variant_id'] = sig.chr.astype(str) + '_' + sig.pos.astype(str) + '_' + sig.ref + '_' + sig.alt

remaining_variants = sig.sort_values('pval').copy()
finemap_variant = [] # lead SNPs to finemap regions around

iteration = 0
while len(remaining_variants) > 0:
    iteration += 1
    logging.info('Iteration {}; {:,} variants remain.'.format(iteration, len(remaining_variants)))
    next_addition = remaining_variants.variant_id.values[0]
    chrom = remaining_variants.chr.values[0]
    pos = remaining_variants.pos.values[0]
    start = pos - EXTENSION
    end = pos + EXTENSION
    finemap_variant.append(next_addition)
    mask = (remaining_variants.chr==chrom) & (remaining_variants.pos>=start) & (remaining_variants.pos<=end)
    remaining_variants = remaining_variants[~mask]

# output regions to finemap
if len(finemap_variant) == 0:
    finemap_regions = pd.DataFrame(columns=['chrom', 'start', 'end', 'variant_id'])
    finemap_regions.to_csv(sys.stdout, sep='\t', index=False, header=False)
else:
    finemap_regions = pd.Series(finemap_variant).str.split('_', expand=True)[[0, 1]]
    finemap_regions.columns = ['chrom', 'pos']
    finemap_regions['end'] = finemap_regions.pos.astype(int)
    finemap_regions['start'] = finemap_regions.end - 1
    finemap_regions['variant_id'] = finemap_variant
    finemap_regions = finemap_regions[['chrom', 'start', 'end', 'variant_id']]
    finemap_regions.start = finemap_regions.start - EXTENSION
    finemap_regions.start = np.maximum(finemap_regions.start, 0)
    finemap_regions.end = finemap_regions.end + EXTENSION
    finemap_regions.to_csv(sys.stdout, sep='\t', index=False, header=False)