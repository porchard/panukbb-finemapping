ROOT=/net/topmed11/working/porchard/panukbb-finemapping
DATA=$(ROOT)/data
WORK=$(ROOT)/work
BIN=$(ROOT)/bin

ANALYSIS=$(WORK)/$@

.PHONY: all

define NL


endef

##### DATA #####
data: variants manifest chain fasta chrom-sizes

variants:
	mkdir -p $(DATA)/$@
	cd $(DATA)/$@ && wget https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_release/full_variant_qc_metrics.txt.bgz

chain:
	mkdir -p $(DATA)/$@/
	cd $(DATA)/$@ && wget http://hgdownload.cse.ucsc.edu/goldenpath/hg38/liftOver/hg38ToHg19.over.chain.gz
	cd $(DATA)/$@ && wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz

chrom-sizes:
	mkdir -p $(DATA)/$@
	cd $(DATA)/$@ && wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes && cat hg19.chrom.sizes | sort -k1,1 > hg19.chrom_sizes && rm hg19.chrom.sizes

manifest:
	mkdir -p $(DATA)/$@
	cd $(DATA)/$@ && wget https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_release/phenotype_manifest.tsv.bgz && zcat phenotype_manifest.tsv.bgz > manifest.tsv && rm phenotype_manifest.tsv.bgz

fasta:
	mkdir -p $(DATA)/$@
	

##### ANALYSES #####

# PanUKBB GWAS fine-mapping
ld-matrices: ANALYSIS=$(WORK)/ld-matrices
ld-matrices:
	mkdir -p $(ANALYSIS)
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results $(ROOT)/download-panukbb-ld-matrices.nf &

lift-variants: ANALYSIS=$(WORK)/lift/variants
lift-variants:
	mkdir -p $(ANALYSIS)
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --variant_file $(DATA)/variants/full_variant_qc_metrics.txt.bgz --chain_file $(DATA)/chain/hg19ToHg38.over.chain.gz --fasta $(DATA)/fasta/Homo_sapiens_assembly38_noALT_noHLA_noDecoy_ERCC.fasta $(ROOT)/lift-panukbb-variants.nf &

get-variants-that-lift-to-same-position:
	mkdir -p $(ANALYSIS)
	zcat $(WORK)/lift/variants/results/vcf/panukbb-hg38.vcf.gz | grep -v "^#" | cut -f1,2,4,5 | sort | uniq -d > $(ANALYSIS)/dups.txt
	cat $(ANALYSIS)/dups.txt | perl -pe 's/\t/_/g' > $(ANALYSIS)/dup-variants.txt
	# grep -f $(ANALYSIS)/dup-variants.txt $(WORK)/ancestry-specific-finemapping/lift-susie/results/susie-cs-and-convergence/gwas.cs.txt # no hits -- so none of the dups are pressent in susie CS lifted to hg38
	zcat $(WORK)/lift/variants/results/conversion-table/snp-names.txt.gz | cut -f4,7 | grep -f $(ANALYSIS)/dup-variants.txt > $(ANALYSIS)/dup-variants.hg19.txt
	cat $(ANALYSIS)/dup-variants.hg19.txt | cut -f1 | perl -pe 's/^/chr/; s/:/_/' > $(ANALYSIS)/dup-variants.hg19.reformatted.txt

select-traits-per-ancestry: ANALYSIS=$(WORK)/selected-traits
select-traits-per-ancestry: ANCESTRIES=AFR EUR
select-traits-per-ancestry:
	mkdir -p $(ANALYSIS)
	$(foreach a,$(ANCESTRIES),python $(BIN)/select-panukbb-gwas-traits-for-ancestry.py $(DATA)/manifest/manifest.tsv $(a) > $(ANALYSIS)/manifest.$(a).txt$(NL))

gwas: ANALYSIS=$(WORK)/gwas
gwas:
	mkdir -p $(ANALYSIS)
	rm -rf $(ANALYSIS)/download.txt
	$(foreach f,$(shell ls $(WORK)/selected-traits/manifest.*.txt),$(BIN)/cut-name aws_path $(f) >> $(ANALYSIS)/download.txt$(NL))
	$(foreach f,$(shell ls $(WORK)/selected-traits/manifest.*.txt),$(BIN)/cut-name aws_path_tabix $(f) >> $(ANALYSIS)/download.txt$(NL))
	cat $(ANALYSIS)/download.txt | sort | uniq | grep -v -w aws_path | grep -v -w aws_path_tabix > $(ANALYSIS)/download.tmp && mv $(ANALYSIS)/download.tmp $(ANALYSIS)/download.txt
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --url_list $(ANALYSIS)/download.txt $(ROOT)/download-panukbb-gwas.nf &

# TODO: update to handle the rare cases where two hg19 variants lift to the same hg38 position
lift-gwas: ANALYSIS=$(WORK)/lift/gwas
lift-gwas:
	mkdir -p $(ANALYSIS)
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --gwas_glob '$(WORK)/gwas/results/*.tsv.bgz' --snp_conversions '$(WORK)/lift/variants/results/conversion-table/snp-names.txt.gz' $(ROOT)/lift-panukbb-gwas.nf &

find-regions-to-finemap: ANALYSIS=$(WORK)/ancestry-specific-finemapping/find-regions-to-finemap
find-regions-to-finemap: ANCESTRIES=AFR EUR
find-regions-to-finemap:
	mkdir -p $(ANALYSIS)
	rm -rf $(ANALYSIS)/traits.txt
	$(foreach a,$(ANCESTRIES),$(BIN)/cut-name trait_id $(WORK)/selected-traits/manifest.$(a).txt | perl -pe 's/$$/\t$(a)/' | awk 'NR>1' >> $(ANALYSIS)/traits.txt$(NL))
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --chrom_sizes $(DATA)/chrom-sizes/hg19.chrom_sizes --gwas_glob '$(WORK)/gwas/results/*.tsv.bgz' --traits $(ANALYSIS)/traits.txt $(ROOT)/panukbb-find-regions-to-finemap.nf &

subset-ld-matrices: ANALYSIS=$(WORK)/ancestry-specific-finemapping/ld-matrices
subset-ld-matrices: ANCESTRIES=AFR EUR
subset-ld-matrices:
	mkdir -p $(ANALYSIS)
	rm -rf $(ANALYSIS)/regions.bed
	$(foreach a,$(ANCESTRIES),cat $(WORK)/ancestry-specific-finemapping/find-regions-to-finemap/results/regions-to-finemap-size-filter/*___$(a).regions.txt | perl -pe 's/$$/\t$(a)/' | sort | uniq | sort -k1,1 -k2n,2 >> $(ANALYSIS)/regions.bed$(NL))
	cat $(ANALYSIS)/regions.bed | sort -k1,1 -k2n,2 > $(ANALYSIS)/regions.tmp && mv $(ANALYSIS)/regions.tmp $(ANALYSIS)/regions.bed
	cd $(ANALYSIS) && nohup nextflow run -resume --pyspark_config '$(ROOT)/config-pyspark.sh' --results $(ANALYSIS)/results --ld_mat_glob '$(WORK)/ld-matrices/results/full-matrix/*.bm' --ld_mat_idx_glob '$(WORK)/ld-matrices/results/full-index/*.ht' --regions $(ANALYSIS)/regions.bed $(ROOT)/panukbb-subset-ld-matrices.nf &

susie: ANALYSIS=$(WORK)/ancestry-specific-finemapping/susie
susie:
	mkdir -p $(ANALYSIS)
	zcat $(WORK)/lift/variants/results/vcf/rejected.vcf.gz | grep -v "^#" | cut -f1,2,4,5 | perl -pe 's/\t/_/g; s/chr//' > $(ANALYSIS)/exclude-snps.txt # variants that fail to lift hg19 --> hg38
	cat $(WORK)/get-variants-that-lift-to-same-position/dup-variants.hg19.reformatted.txt | perl -pe 's/^chr//' >> $(ANALYSIS)/exclude-snps.txt
	cat $(ANALYSIS)/exclude-snps.txt | sort > $(ANALYSIS)/exclude-snps.tmp && mv $(ANALYSIS)/exclude-snps.tmp $(ANALYSIS)/exclude-snps.txt
	python $(BIN)/make-ancestry-specific-finemapping-susie-input-list.py --gwas-glob '$(WORK)/gwas/results/*.bgz' --regions-glob '$(WORK)/ancestry-specific-finemapping/find-regions-to-finemap/results/regions-to-finemap-size-filter/*' --ld-matrix-glob '$(WORK)/ancestry-specific-finemapping/ld-matrices/results/ld-matrix/*.ld.tsv.bgz' --ld-matrix-variants-glob '$(WORK)/ancestry-specific-finemapping/ld-matrices/results/ld-matrix/*.variants.tsv.gz' | awk 'NR>1' > $(ANALYSIS)/to-run.txt
	cd $(ANALYSIS) && nohup nextflow run -resume -with-trace --to_run $(ANALYSIS)/to-run.txt --snps_exclude $(ANALYSIS)/exclude-snps.txt --results $(ANALYSIS)/results $(ROOT)/run-susieR.nf &

filter-susie: ANALYSIS=$(WORK)/ancestry-specific-finemapping/filter-susie
filter-susie:
	mkdir -p $(ANALYSIS)/susie-out
	cd $(ANALYSIS) && python $(BIN)/remove-regions-with-extreme-number-of-variants.py $(WORK)/ancestry-specific-finemapping/susie/to-run.txt '$(WORK)/ancestry-specific-finemapping/ld-matrices/results/ld-matrix/*.variants.tsv.gz' $(WORK)/ancestry-specific-finemapping/susie/results/susie-out $(ANALYSIS)/susie-out

lift-susie: ANALYSIS=$(WORK)/ancestry-specific-finemapping/lift-susie
lift-susie:
	mkdir -p $(ANALYSIS)/data
	cd $(ANALYSIS) && nohup nextflow run -resume --susie_glob '$(WORK)/ancestry-specific-finemapping/filter-susie/susie-out/*' --to_run $(WORK)/ancestry-specific-finemapping/filter-susie/to-run.txt --chain $(DATA)/chain/hg19ToHg38.over.chain.gz --snp_conversions $(WORK)/lift/variants/results/conversion-table/snp-names.txt.gz --snp_conversions_idx $(WORK)/lift/variants/results/conversion-table/snp-names.txt.gz.tbi --results $(ANALYSIS)/results $(ROOT)/lift-susieR.nf &