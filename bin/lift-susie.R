#!/usr/bin/env Rscript

library(susieR)
library(dplyr)
library(glue)


tabix <- function(f, region) {
       #region <- 'chr1:149964027-150464028'
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

#args <- c('snp-names.txt.gz', 'hg38-regions.bed', 'susie-list.txt')
args <- commandArgs(T)
SNP_CONVERSIONS_FILE <- args[1]
LIFTOVER_FILE <- args[2]
SUSIE_LIST <- args[3]


susie_list <- read.table(SUSIE_LIST, head=F, as.is=T, col.names=c('rda'))
susie_list$hg19_region <- sapply(strsplit(susie_list$rda, '___'), function(x){gsub('.rda', '', x[3])})
susie_list$hg19_region <- paste0('chr', gsub('_', ':', susie_list$hg19_region))

liftover <- read.table(LIFTOVER_FILE, head=F, as.is=T, col.names=c('chrom', 'start', 'end', 'hg19'))
liftover$hg38 <- with(liftover, paste(chrom, start, end, sep=':'))
region_conversions <- liftover$hg38
names(region_conversions) <- liftover$hg19
susie_list$hg38_region <- region_conversions[susie_list$hg19_region]

susie_list$rda_hg38 <- NA
susie_list$REASON_FAIL <- NA
susie_list$TOTAL_SNPS <- NA
susie_list$SNPS_LIFTED <- NA

for(i in 1:nrow(susie_list)) {
    hg38_region <- susie_list$hg38_region[i]
    if (is.na(hg38_region)) {
        susie_list$REASON_FAIL[i] <- 'region_liftover_failed'
        next
    }
    conversions <- tabix(SNP_CONVERSIONS_FILE, sub('-', ':', gsub(':', '-', hg38_region)))

    hg19_to_hg38 <- conversions$hg38_id
    names(hg19_to_hg38) <- gsub('chr', '', conversions$hg19_id)

    rda_in <- susie_list$rda[i]
    rda_out <- paste(c(strsplit(rda_in, '___')[[1]][1:2], hg38_region), collapse='___')
    rda_out <- paste0(rda_out, '.rda')
    rda_out <- gsub(':', '_', rda_out)
    rda_out <- basename(rda_out)
    susie_list$rda_hg38[i] <- rda_out

    load(rda_in) # named susie_out

    # missing some conversions?
    missing <- !(colnames(susie_out$alpha) %in% names(hg19_to_hg38))

    susie_list$TOTAL_SNPS[i] <- length(missing)
    susie_list$SNPS_LIFTED[i] <- sum(!missing)

    index_and_snp_conversions <- data.frame(
    hg19_snp=colnames(susie_out$alpha),
    hg19_index=seq(1, length(colnames(susie_out$alpha))),
    hg38_snp=as.vector(hg19_to_hg38[colnames(susie_out$alpha)])
    )
    index_and_snp_conversions$hg38_index <- NA
    index_and_snp_conversions$hg38_index[!is.na(index_and_snp_conversions$hg38_snp)] <- seq(1, sum(!is.na(index_and_snp_conversions$hg38_snp)))

    hg19_to_hg38_index <- index_and_snp_conversions$hg38_index
    names(hg19_to_hg38_index) <- index_and_snp_conversions$hg19_index

    susie_out$alpha <- susie_out$alpha[,colnames(susie_out$alpha) %in% names(hg19_to_hg38)] # remove missing elements
    colnames(susie_out$alpha) <- hg19_to_hg38[colnames(susie_out$alpha)] # convert

    susie_out$mu <- susie_out$mu[,colnames(susie_out$mu) %in% names(hg19_to_hg38)]
    colnames(susie_out$mu) <- hg19_to_hg38[colnames(susie_out$mu)]

    # change mu sign if the ref/alt alleles flip between refs
    # doing this after converting to hg38 because otherwise some may be missing
    swapped_alleles <- (conversions$swapped_alleles=='True')
    names(swapped_alleles) <- conversions$hg38_id
    mu_sign_factor <- ifelse(swapped_alleles[colnames(susie_out$mu)] == F, 1, -1)
    stopifnot(names(mu_sign_factor) == colnames(susie_out$mu))
    if (any(mu_sign_factor == -1)) {
        new_mu <- t(t(susie_out$mu) * mu_sign_factor)
    }


    susie_out$mu2 <- susie_out$mu2[,colnames(susie_out$mu2) %in% names(hg19_to_hg38)]
    colnames(susie_out$mu2) <- hg19_to_hg38[colnames(susie_out$mu2)]

    susie_out$lbf_variable <- susie_out$lbf_variable[,colnames(susie_out$lbf_variable) %in% names(hg19_to_hg38)]
    colnames(susie_out$lbf_variable) <- hg19_to_hg38[colnames(susie_out$lbf_variable)]

    susie_out$pi <- susie_out$pi[!missing]
    susie_out$X_column_scale_factors <- susie_out$X_column_scale_factors[!missing]

    susie_out$pip <- susie_out$pip[!missing]
    names(susie_out$pip) <- hg19_to_hg38[names(susie_out$pip)]

    susie_out$Rr <- NULL

    sets_new_cs <- lapply(susie_out$sets$cs, function(x){
        y <- hg19_to_hg38_index[x]
        y <- y[!is.na(y)]
        return(y)
    })

    susie_out$sets$cs <- sets_new_cs
    save(susie_out, file=rda_out)

}


write.table(susie_list, file='susie-conversion-summary.txt', sep='\t', quote=F, append=F, row.names=F)