#!/usr/bin/env python
# coding: utf-8

import sys
import glob
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def symlink(input, output):
    if not os.path.exists(output):
        os.symlink(input, output)
    else:
        raise ValueError('Output file already exists: {}'.format(output))


TO_RUN, VARIANT_LIST_GLOB, INDIR, OUTDIR = sys.argv[1:]
#TO_RUN = '/net/topmed11/working/porchard/panukbb-finemapping/work/ancestry-specific-finemapping/susie/to-run.txt'
#INDIR = '/net/topmed11/working/porchard/panukbb-finemapping/work/ancestry-specific-finemapping/susie/results/susie-out'
#VARIANT_LIST_GLOB = '/net/topmed11/working/porchard/panukbb-finemapping/work/ancestry-specific-finemapping/ld-matrices/results/ld-matrix/*.variants.tsv.gz'
#OUTDIR = '/net/topmed11/working/porchard/panukbb-finemapping/work/filter-susie/susie-out'
THRESHOLD = 40000


to_run = pd.read_csv(TO_RUN, sep='\t', header=None, names=['trait', 'gwas', 'gwas_index', 'ancestry', 'chrom', 'start', 'end', 'ld', 'variants'])
to_run['output'] = to_run.trait + '___' + to_run.ancestry + '___' + to_run.chrom.astype(str) + '_' + to_run.start.astype(str) + '_' + to_run.end.astype(str) + '.rda'

# determine number of SNPs in each region to finemap
variant_lists = glob.glob(VARIANT_LIST_GLOB)

n_variants = [[os.path.basename(f), len(pd.read_csv(f, sep='\t'))] for f in variant_lists]
n_variants = pd.DataFrame(n_variants, columns=['file', 'n_variants'])
n_variants[['chrom', 'start', 'end', 'ancestry']] = n_variants.file.str.extract('(.*)_(.*)_(.*)___(.*)\.variants\.tsv\.gz', expand=True)
n_variants.chrom = n_variants.chrom.astype(str)
n_variants.start = n_variants.start.astype(int)
n_variants.end = n_variants.end.astype(int)
n_variants = n_variants[['chrom', 'start', 'end', 'ancestry', 'n_variants']]

df = to_run.merge(n_variants)
df['region_length'] = df.end - df.start


fig, axs = plt.subplots(ncols=2, figsize=(10, 5))

ax = axs[0]
df.n_variants.hist(bins=100, ax=ax)
ax.set_xlabel('Number of variants in region')
ax.set_ylabel('Number of regions')
ax.axvline(THRESHOLD, color='red', ls='--')
ax.set_title('{:,} of {:,} regions have > {:,} variants'.format((df.n_variants > THRESHOLD).sum(), df.shape[0], THRESHOLD))

ax = axs[1]
sns.scatterplot(y='n_variants', x='region_length', data=df, ax=ax, alpha=0.1)
ax.set_xlabel('Region length')
ax.set_ylabel('Number of variants in region')
ax.axhline(THRESHOLD, color='red', ls='--')

fig.tight_layout()
fig.savefig('summary.png', dpi=300, bbox_inches='tight', facecolor='white')
fig.clf()

df.loc[df.n_variants <= THRESHOLD,to_run.columns.to_list()].drop(columns=['output']).to_csv('to-run.txt', sep='\t', index=False, header=False)

# symlink the regions with <= THREHSOLD variants
to_symlink = df[df.n_variants <= THRESHOLD]

for f in to_symlink.output:
    input = os.path.join(INDIR, f)
    output = os.path.join(OUTDIR, f)
    symlink(input, output)