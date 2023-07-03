#!/usr/bin/env Rscript
# coding: utf-8

library(susieR)
library(glue)

args <- commandArgs(T)

GWAS <- args[1]
ANCESTRY <- args[2]
LD <- args[3]
LD_VARIANTS <- args[4]
REGION <- args[5] # chr:1-1000
IGNORE_SNPS <- args[6]
OUT <- args[7]
DIAGNOSTIC <- args[8]

tabix <- function(f, region) {
    # region follows format '{chrom}:{start}-{end}'
    header <- as.character(as.vector(read.table(gzfile(f), nrows=1, sep='\t')[c(1),]))
    x <- strsplit(system(glue('tabix {f} {region}'), intern = TRUE,
        ignore.stdout = FALSE, ignore.stderr = FALSE,
        wait = TRUE, timeout = 0), '\t')
    y <- list()
    for(i in 1:length(header)) {
        y[[i]] <- sapply(x, function(z){z[i]})
    }
    names(y) <- header
    df <- data.frame(y)
    return(df)
}


EXCLUDE_SNPS <- read.table(IGNORE_SNPS, sep='\t', head=F, as.is=T)$V1

# load in the GWAS

gwas <- tabix(GWAS, REGION)
gwas$variant_id <- with(gwas, paste(chr, pos, ref, alt, sep='_'))
gwas <- gwas[,c('variant_id', glue('beta_{ANCESTRY}'), glue('se_{ANCESTRY}'), glue('low_confidence_{ANCESTRY}'))]
colnames(gwas) <- c('variant_id', 'beta', 'se', 'low_confidence')
gwas <- gwas[gwas$low_confidence=='false',]
gwas$beta <- as.numeric(gwas$beta)
gwas$se <- as.numeric(gwas$se)
gwas$z <- gwas$beta / gwas$se
gwas <- gwas[!is.na(gwas$beta),]
NUMBER_SNPS_EXCLUDE <- sum(gwas$variant_id %in% EXCLUDE_SNPS)
if (NUMBER_SNPS_EXCLUDE > 0) {
    print(glue("Excluding {NUMBER_SNPS_EXCLUDE} variants from {IGNORE_SNPS}"))
    gwas <- gwas[!gwas$variant_id %in% EXCLUDE_SNPS,]
}


ld <- read.table(gzfile(LD), sep='\t', head=F, colClasses='numeric', row.names=NULL)
ld <- as.matrix(ld)
ld_variants <- read.table(gzfile(LD_VARIANTS), sep='\t', head=T)
stopifnot(length(ld_variants$variant_id) == nrow(ld))
rownames(ld) <- ld_variants$variant_id
colnames(ld) <- ld_variants$variant_id

isnps <- intersect(gwas$variant_id, ld_variants$variant_id)
NDROPPED <- nrow(gwas) - length(isnps)
NTOTAL <- nrow(gwas)
warning(glue('Dropping {NDROPPED} of {NTOTAL} GWAS SNPs (missing from LD matrix)'))
gwas <- gwas[gwas$variant_id %in% isnps,]

ld <- ld[gwas$variant_id,gwas$variant_id]
# currently upper triangular matrix; make symmetric
symm <- ld + t(ld)
diag(symm) <- diag(symm) / 2


n_cs <- function(susie_obj) {
    return(length(susie_obj$sets$cs))
}


run_susie_and_select_l <- function(z, R, variant_ids, L_consider=c(10, 20, 30, 40)) {
    L_SORTED <- sort(unique(L_consider))
    N_CS <- c()
    SUSIE_OUT <- list()
    # run with smallest L. if # CS discovered >= 0.7*L, try next L. Repeat until # CS discovered < 0.7*L. Then select min(L[L>max_cs_discovered])
    for(L in L_SORTED) {
        print(glue('Running w/ L = {L}'))
        susie_out_l <- susie_rss(z = z, R = symm, L = L)
        colnames(susie_out_l$lbf_variable) <- variant_ids
        colnames(susie_out_l$alpha) <- variant_ids
        names(susie_out_l$pip) <- variant_ids
        SUSIE_OUT[[length(SUSIE_OUT)+1]] <- susie_out_l
        if (!susie_out_l$converged){
            print('Failed to converge; trying next L')
            N_CS <- c(N_CS, 0)
            next
        }
        cs_l <- n_cs(susie_out_l)
        print(glue('Discovered {cs_l} CS'))
        N_CS <- c(N_CS, cs_l)
        if (cs_l < 0.7*L) {
            break
        }
    }
    max_cs_discovered <- max(N_CS)
    L <- min(L_SORTED[L_SORTED>max_cs_discovered])
    IDX <- which(L_SORTED==L)
    print(glue('Selected L = {L}'))
    return(SUSIE_OUT[[IDX]])
}

# plot diagnostic
png(DIAGNOSTIC, height=4, width=4, units="in", res=150)
condz_in = kriging_rss(gwas$z, symm)
print(condz_in$plot)
dev.off()

susie_out <- run_susie_and_select_l(gwas$z, symm, gwas$variant_id, c(10, 20, 30, 40))
save(susie_out, file=OUT)
