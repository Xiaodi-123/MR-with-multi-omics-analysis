load('table1.Rdata')
devtools::install_github("phenoscanner/phenoscanner")
install.packages('MendelianRandomization')
library (MendelianRandomization)
installed.packages('phenoscanner')
BiocManager::install('phenoscanner')
library(phenoscanner)
mr_phenoscanner <- function(dat, 
                            catalog = c("GWAS", "eQTL", "pQTL", "mQTL", "methQTL"),
                            pvalue = 5e-08,
                            proxies = "EUR", 
                            r2 = 0.8,
                            build = 37)
{
  stopifnot("SNP" %in% names(dat))
  
  snpquery <- unique(dat$SNP)
  query_limit <- ceiling(0.01*length(snpquery))
  
  query_res <- list()
  for (i in 1:query_limit) {
    lo_i <- 100*i - 99
    up_i <- min(100*i,length(snpquery))
    snp_i <- snpquery[lo_i:up_i]
    print(paste0("Searching for ", snp_i))
    
    for (catalogue in catalog) {
      query_i <- try(
        phenoscanner::phenoscanner(
          snpquery = snp_i, catalogue = catalogue,
          proxies = proxies, pvalue = pvalue, r2 = r2, build = build)
      )
      if(!"try-error" %in% class(query_i))
      {
        if("results" %in% names(query_i))
        {
          query_i <- list(query_i$results)
          names(query_i) <- catalogue
          query_res <- append(query_res, query_i)
        }
      }
    }
  }
  
  return(query_res)
}


pleio_pheno <- mr_phenoscanner(table1) # the key information is the SNP column


pleio_pheno_summary_raw <- pleio_pheno$GWAS %>% 
  mutate(catalog = "GWAS") %>% 
  filter(p != "NA",
         ancestry %in% c("Mixed", "European")) %>% 
  mutate(across(p, as.numeric)) %>% 
  group_by(trait, pmid) %>% 
  filter(p == min(p)) %>% 
  ungroup()


