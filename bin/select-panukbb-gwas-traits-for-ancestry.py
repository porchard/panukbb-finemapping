#!/usr/bin/env python
# coding: utf-8


import sys
import pandas as pd


MANIFEST_FULL, ANCESTRY_OF_INTEREST = sys.argv[1:]
ANCESTRIES = ['AFR', 'AMR', 'CSA', 'EAS', 'EUR', 'MID']


manifest = pd.read_csv(MANIFEST_FULL, sep='\t')

# get columns referring to each ancestry
ancestry_columns = {ancestry: manifest.columns.to_series().str.contains(f'_{ancestry}$', regex=True).where(lambda x: x==True).dropna().index.to_list() for ancestry in ANCESTRIES}

# drop columns referring to ancestries that aren't our ancestry of interest
for ancestry, columns in ancestry_columns.items():
    if ancestry != ANCESTRY_OF_INTEREST:
        manifest = manifest.drop(columns=columns)


manifest['trait_id'] = manifest.filename.str.replace('.tsv.bgz', '', regex=False)
manifest.description = manifest.description.fillna('')
manifest.category = manifest.category.fillna('')

# keep pass only, and max independent set
manifest = manifest[manifest[f'phenotype_qc_{ANCESTRY_OF_INTEREST}']=='PASS']
manifest = manifest[manifest.in_max_independent_set]

# drop other unnecessary columns
drop_columns = ['filename', 'filename_tabix', 'pops', 'num_pops', 'num_pops_pass_qc', 'pops_pass_qc']
drop_columns += [i for i in manifest.columns if 'n_cases_hq_cohort' in i]
drop_columns += [i for i in manifest.columns if 'n_cases_full_cohort' in i]
manifest = manifest.drop(columns=drop_columns)

# remove any columns w/ only a single value
manifest = manifest.loc[:,(manifest.nunique()>1).to_list()]


# drop certain selected traits
DROP_DESCRIPTION = [
    'reproduciblity of spirometry measurement using ERS/ATS criteria',
    'leisure/social activities',
    'recent feelings'
]

DROP_CATEGORY = [
    'mental distress',
    'mental health',
    'mental disorders',
    'sexual',
    'education',
    'intelligence',
    'alcohol'
]


DROPPED = []
for i in DROP_DESCRIPTION:
    DROP = manifest[manifest.description.str.contains(i, case=False, regex=False)]
    DROPPED.append(DROP)
    manifest = manifest[~manifest.description.str.contains(i, case=False, regex=False)]
for i in DROP_CATEGORY:
    DROP = manifest[manifest.category.str.contains(i, case=False, regex=False)]
    DROPPED.append(DROP)
    manifest = manifest[~manifest.category.str.contains(i, case=False, regex=False)]
DROPPED = pd.concat(DROPPED)

# remove pipes in some of the names
manifest.phenocode = manifest.phenocode.str.replace('|', '_', regex=False)
manifest.trait_id = manifest.trait_id.str.replace('|', '_', regex=False)

manifest.to_csv(sys.stdout, sep='\t', index=False)