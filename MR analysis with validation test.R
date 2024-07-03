library(tidyverse)

load('exposure_dat.Rdata')

gene=read.table('./OR.txt',sep = '\t',header = T)
gene=gene$exposure

exposure_dat=exposure_dat[exposure_dat$exposure %in% gene,]
library(TwoSampleMR)
outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes="ebi-a-GCST90018863")
exposure_dat=exposure_dat[exposure_dat$SNP %in% outcome_dat$SNP,]
harmonised_dat <- harmonise_data(exposure_dat, outcome_dat)
mr_modified <- function(dat = harmonised_dat, prop_var_explained = T)
{
  mr_res <- mr(dat)
  
  pve <- dat %>% 
    dplyr::select(id.exposure, beta.exposure, se.exposure, samplesize.exposure) %>% 
    dplyr::group_by(id.exposure) %>% 
    dplyr::summarise(pve = sum((beta.exposure^2)/(beta.exposure^2 + samplesize.exposure*se.exposure^2)))
  
  if(prop_var_explained)
  {
    mr_res <- mr_res %>% 
      dplyr::left_join(pve, by = "id.exposure")
  }
  
  return(mr_res)
}
mr_res_vali <- mr_modified(harmonised_dat, prop_var_explained = T)

save(mr_res_vali,harmonised_dat,file ='mr_input_res_vali.Rdata')


result_or=generate_odds_ratios(mr_res_vali)
write.table(result_or[,4:ncol(result_or)],"OR_vali.txt",row.names = F,sep = "\t",quote = F)
library(grid)
library(forestploter)
mydata=read.table("OR_vali.txt",header = T,sep = "\t")
mydata$` ` <- paste(rep(" ", 20), collapse = " ")
mydata$`OR (95% CI)` <- ifelse(is.na(mydata$or), "",sprintf("%.4f (%.4f - %.4f)",
                                                            mydata$or, mydata$or_lci95, 
                                                            mydata$or_uci95))
forest(mydata[,c(1:3,6,13,14)],
       est = mydata$or,
       lower =mydata$or_lci95, 
       upper = mydata$or_uci95,
       sizes =0.3,
       ci_column =5 ,
       ref_line = 1,
       xlim = c(0.05, 3),
)
mr_res1=mr_res_vali[mr_res_vali$exposure=='N4BP2L1',]
harmonised_dat1=harmonised_dat[harmonised_dat$exposure=='N4BP2L1',]

mr_scatter_plot(mr_res1, harmonised_dat1)

res_single=mr_singlesnp(harmonised_dat1)
mr_forest_plot(res_single)

mr_funnel_plot(singlesnp_results = res_single)

mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(harmonised_dat1))
