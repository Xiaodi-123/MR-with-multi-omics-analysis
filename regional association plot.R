devtools::install_github("boxiangliu/locuscomparer")
library(locuscomparer)

snp_overlap <- intersect(my_eqtl[['snp']],
                         gwas[['snp']])
snp_overlap <- unique(snp_overlap)

exposure_dat <-my_eqtl[my_eqtl[['snp']] %in% snp_overlap,]
outcome_dat <- gwas[gwas[['snp']] %in% snp_overlap,]

exposure_dat <- exposure_dat[order(exposure_dat[['snp']]),]
outcome_dat <- outcome_dat[order(outcome_dat[['snp']]),]

exposure_dat <- data.frame(
  rsid = exposure_dat[['snp']],
  pval = exposure_dat[['pvalues']])

outcome_dat <- data.frame(
  rsid = outcome_dat[['snp']],
  pval = outcome_dat[['pvalues']])



write.table(exposure_dat,"coloc_exposure_test.tsv",sep = "\t",row.names = F,quote = F)
write.table(outcome_dat,"coloc_outcome_test.tsv",sep = "\t",row.names = F,quote = F)

library(locuscomparer)
p <- locuscompare(in_fn1 ="coloc_outcome_test.tsv", 
                  in_fn2 ="coloc_exposure_test.tsv", 
                  title1 = paste0('ILD'," GWAS"),  
                  title2 = paste0('N4BP2L1'," eQTL"), 
                  combine =T
)
p

