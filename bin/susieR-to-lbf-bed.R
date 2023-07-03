#!/usr/bin/env Rscript

library(susieR)
library(dplyr)

args <- commandArgs(T)
TRAIT <- args[1]
OUT <- args[2]
INPUT <- args[3:length(args)]

results <- list()

for(i in INPUT) {
    tmp <- load(i)
    s <- get(tmp)
    if (!s$converged) {
        next
    }

    lbf <- s$lbf_variable
    rownames(lbf) <- paste0('L', 1:nrow(lbf))
    lbf <- as.data.frame(t(lbf))
    EFFECTS <- colnames(lbf)

    variant_info <- strsplit(rownames(lbf), '_')
    chrom <- sapply(variant_info, function(x){x[1]})
    pos <- as.numeric(sapply(variant_info, function(x){x[2]}))
    lbf$chrom <- chrom
    lbf$start <- pos - 1
    lbf$end <- pos
    lbf$variant_id <- rownames(lbf)
    lbf$trait <- TRAIT
    NONEFFECT_COLS <- colnames(lbf)[!colnames(lbf) %in% EFFECTS]
    lbf <- lbf[,c(NONEFFECT_COLS, EFFECTS)]
    colnames(lbf)[colnames(lbf)=='chrom'] <- '#chrom'
    results[[length(results)+1]] <- lbf

}

results <- bind_rows(results)
write.table(results, OUT, sep='\t', row.names=F, quote=F)
