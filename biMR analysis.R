bimr_ILD<- extract_instruments('ebi-a-GCST90018643', p1=5e-08, clump=TRUE)

outcome_gene<- extract_outcome_data(snps=bimr_ILD$SNP, outcomes="eqtl-a-ENSG00000139597")

bimr_ILD =bimr_ILD [bimr_ILD $SNP %in% outcome_gene$SNP,]

harmonised_ILD_gene <- harmonise_data(bimr_ILD, outcome_gene)

bimr_mr_ILD_gene <- mr(harmonised_ILD_gene)

result_or=generate_odds_ratios(bimr_mr_ILD_gene)
write.table(result_or[,4:ncol(result_or)],"bi_OR.txt",row.names = F,sep = "\t",quote = F)

library(grid)
library(forestploter)

mydata=read.table("bi_OR.txt",header = T,sep = "\t")
mydata$outcome='N4BP2L1'
mydata$` ` <- paste(rep(" ", 20), collapse = " ")
mydata$`OR (95% CI)` <- ifelse(is.na(mydata$or), "",sprintf("%.4f (%.4f - %.4f)",
                                                            mydata$or, mydata$or_lci95, 
                                                            mydata$or_uci95))
forest(mydata[,c(1:3,6,12,13,14)],
       est = mydata$or,
       lower =mydata$or_lci95, 
       upper = mydata$or_uci95,
       sizes =0.3,
       ci_column =6 ,
       ref_line = 1,
       xlim = c(0.05, 3),
)