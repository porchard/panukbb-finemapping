#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sys
import logging

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s: %(message)s')
VARIANT_LIFTOVER_INFO = sys.argv[1]

COMPLEMENTS = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

def reverse_complement(x):
    return ''.join([COMPLEMENTS[i] for i in x[::1]])


# In[ ]:


print('\t'.join(['chrom', 'start', 'end', 'rsid', 'ref', 'alt', 'hg38_id', 'hg19_id', 'reverse_complemented_alelles', 'swapped_alleles']))
#skip_count = 0
#total_count = 0
with open(VARIANT_LIFTOVER_INFO, 'r') as fh:
    for line in fh:
        #total_count += 1
        ALLELE_MISMATCH = False
        chrom, pos, rsid, ref, alt, variant_info = line.rstrip().split()
        end = pos
        start = str(int(end) - 1)
        variant_info = variant_info.split(';')
        hg19_chrom = None
        hg19_pos = None
        hg19_ref = ref
        hg19_alt = alt
        reverse_complemented_alelles = 'ReverseComplementedAlleles' in variant_info
        swapped_alleles = 'SwappedAlleles' in variant_info
        if reverse_complemented_alelles:
            hg19_ref = reverse_complement(hg19_ref)
            hg19_alt = reverse_complement(hg19_alt)
        if swapped_alleles:
            (hg19_ref, hg19_alt) = (hg19_alt, hg19_ref)
        for i in variant_info:
            if 'OriginalAlleles' in i:
                original_alleles = i.replace('OriginalAlleles=', '').split(',')
                if hg19_ref not in original_alleles or hg19_alt not in original_alleles:
                    ALLELE_MISMATCH = True
                    #skip_count += 1
                    logging.warning(f'Skipping variant {rsid} ({line.rstrip()}); original alleles dont match inferred original alleles')
                continue
            if 'OriginalContig' in i:
                hg19_chrom = i.replace('OriginalContig=', '')
                continue
            if 'OriginalStart' in i:
                hg19_pos = i.replace('OriginalStart=', '')
                continue
        for i in [hg19_chrom, hg19_pos, hg19_ref, hg19_alt]:
            if i is None:
                raise ValueError(line)
        hg38_id = f'{chrom}_{pos}_{ref}_{alt}'
        hg19_id = f'{hg19_chrom}_{hg19_pos}_{hg19_ref}_{hg19_alt}' if not ALLELE_MISMATCH else '-'
        #if not SKIP:
        print('\t'.join([chrom, start, end, rsid, ref, alt, hg38_id, hg19_id, str(reverse_complemented_alelles), str(swapped_alleles)]))
#logging.info(f'Skipped {skip_count} of {total_count}')
logging.info('Done.')
        
        

