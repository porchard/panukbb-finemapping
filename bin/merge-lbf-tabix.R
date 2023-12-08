#!/usr/bin/env Rscript

library(dplyr)

args <- commandArgs(T)
OUT <- args[1]
IN <- args[2:length(args)]

tmp <- lapply(IN, function(f){
    read.table(f, head=T, sep='\t', as.is=T, comment.char='')
})

combined <- bind_rows(tmp)

colnames(combined) <- gsub('X.', '#', colnames(combined))
options(scipen = 999)
write.table(combined, file=OUT, sep='\t', append=F, quote=F, row.names=F)
