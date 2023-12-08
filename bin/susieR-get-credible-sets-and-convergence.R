#!/usr/bin/env Rscript

library(susieR)
library(dplyr)
library(glue)

args <- commandArgs(T)
SUSIE_LIST <- args[1]
OUT_PREFIX <- args[2]

susie_list <- read.table(SUSIE_LIST, sep='\t', head=F, as.is=T, col.names=c('rda')) # basenames MUST be unique
stopifnot(max(table(basename(susie_list$rda))) == 1) # check that basenames are unique

# get CS locations


all_results <- list()
converged <- c()

for(rda_in in susie_list$rda) {
    tmp <- load(rda_in)
    s <- get(tmp)

    converged <- c(converged, s$converged)
    if (!s$converged) {
        next
    }
    cs <- s$sets$cs
    cs <- cs[lapply(cs, length) > 0]
    if (length(cs) > 0) {
        cs <- lapply(cs, function(x){data.frame(snp=colnames(s$alpha)[x])})
        for(i in 1:length(cs)) {
            cs[[i]]$cs <- names(cs)[i]
            rownames(cs[[i]]) <- paste(cs[[i]]$cs, cs[[i]]$snp)
            # cs[[i]]$pip <- s$pip[cs[[i]]$snp] # use this if want PIP aggregated across effects
            cs[[i]]$pip <- s$alpha[as.numeric(unique(gsub('L', '', names(cs)[i]))),cs[[i]]$snp] # use this if want single effect PIP
        }
        cs <- bind_rows(cs)
        cs$rda <- basename(rda_in)
        cs$converged <- s$converged
        all_results[[length(all_results)+1]] <- cs
    }
    
}

all_results <- bind_rows(all_results)
cs <- all_results[all_results$converged,c('snp', 'cs', 'pip', 'rda')]
converged <- data.frame(rda=basename(susie_list$rda), converged=converged)

write.table(cs, file=glue('{OUT_PREFIX}cs.txt'), sep='\t', col.names=F, row.names=F, quote=F)
write.table(converged, file=glue('{OUT_PREFIX}converged.txt'), sep='\t', col.names=F, row.names=F, quote=F)
