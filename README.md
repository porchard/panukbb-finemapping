# Finemapping of PanUKBB GWAS

## Setup

In `Makefile`, set `ROOT` to the current directory.

Download chain files and some PanUKBB resources:

1. `make data`
2. Fetch the fasta file from [here](https://console.cloud.google.com/storage/browser/_details/gtex-resources/references/Homo_sapiens_assembly38_noALT_noHLA_noDecoy_ERCC.fasta;tab=live_object?pli=1) and place in data/fasta


## Dependencies

You must have the following installed:

1. NextFlow (>= v. 21.04). NextFlow should be configured as appropriate for your computing platform.
2. Singularity (v. 3) installed
3. `aws` CLI
4. python3, with `pandas`, `hail`, `matplotlib`, and `seaborn` packages

## Running

1. Download the LD matrices: `make ld-matrices`
2. Lift PanUKBB variants to hg38: `make lift-variants`
3. Select the traits to be finemapped: `make select-traits-per-ancestry`
4. Download the GWAS: `make gwas`
5. Lift GWAS to hg38: `make lift-gwas`
6. Identify regions to finemap: `make find-regions-to-finemap`
7. Subset LD matrices: `make subset-ld-matrices`
8. Run SuSIE: `make susie`
9. Filter SuSiE output, removing regions w/ more than 40k variants: `make filter-susie`
10. Lift SuSiE to hg38: `make lift-susie`
