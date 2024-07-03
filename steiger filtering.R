load('mr_input_res.Rdata')
steiger_res= harmonised_dat %>%
  filter(SNP %in% table1$SNP) %>% 
  steiger_filtering()
steiger_res$steiger_pval

steiger_res2= harmonised_dat %>%
  filter(SNP %in% table1$SNP) %>% 
  directionality_test()

write.csv(table1,file ='table1.csv',quote = F)

