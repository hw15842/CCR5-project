######################################
### CCR5_PheWAS_for_new_395_traits ###
######################################

setwd("/Users/hw15842/Documents/PhD work/CCR5_paper")

library (devtools)
library(plyr)
library(data.table)

install_github("MRCIEU/TwoSampleMR")

library(TwoSampleMR)

ao <- available_outcomes()

## Exposure CCR5 data ##
exposure_CCR5 <- read.table("CCR5.txt", header=T)

exposure_CCR5 <- format_data(exposure_CCR5, type="exposure")

## Outcome traits that are in MRbase
outcome_traits_mrbase <- read.csv("list_of_395_traits.csv", header=T)


## Outcome traits that are not in MRbase
outcome_traits_other <- fread("other_GWAS_CCR5.txt", header=T)

outcome_traits_other <- format_data(outcome_traits_other, type="outcome")




## Remove NAs from mr.base.ids column in outcome_traits_mrbase ##

outcome_traits_mrbase_na.rm <- outcome_traits_mrbase[!is.na(outcome_traits_mrbase$MR.Base.id),]



## Check which outcomes are giving nulls and remove them from the dataset
outcomes_check_nulls_ <- function(outcome_ID){
  extract_outcome_data(snps = exposure_CCR5$SNP, outcome = outcome_ID )
}

outcomes_with_nulls <- lapply(outcome_traits_mrbase_na.rm$MR.Base.id, outcomes_check_nulls_)
outcomes_with_nulls_2 <- ldply(outcomes_with_nulls, data.frame)

df_new <- sapply(outcomes_with_nulls, is.null)
table(df_new) 
## rows 6 and 9 are null so remove those from the dataset

outcome_traits_mrbase_na.rm <- outcome_traits_mrbase_na.rm[-c(6,9),]

## Run the MRbase outcomes ##

run_mr_CCR5 <- function(outcome_ID){
  
  outcome_1 <- extract_outcome_data(snps = exposure_CCR5$SNP, outcome = outcome_ID)
  dat1 <- harmonise_data(exposure_CCR5, outcome_1)
  res1 <- mr(dat1, method_list = "mr_ivw")
  return(res1)
}

results_MRbase <- lapply(outcome_traits_mrbase_na.rm$MR.Base.id, run_mr_CCR5)
results_MRbase_table <- ldply(results_MRbase, data.frame)

length(results_MRbase)
nrow(results_MRbase_table) ## 4 of them didnt have results 

write.table(results_MRbase_table, file = "results_MRbase.txt", quote=F, sep="\t")


### run the other (non mrbase) traits ##

dat <- harmonise_data(exposure_dat = exposure_CCR5, outcome_dat = outcome_traits_other)

results_other <- mr(dat, method_list=c("mr_ivw"))

write.table(results_other, file="results_other.txt", quote=F, sep="\t")


## Put both sets of results together ##

results_all <- rbind(results_MRbase_table, results_other)

write.table(results_all, file="results_all.txt", quote=F, sep="\t")


## Look for significance ##

results_all_sorted_pval <- results_all[order(results_all$pval),]

0.05/386

## WBC count, RBC count and osteoporosis all significant after adjusting for multiple testing ##


## Pick out WBC and RBC count and osteoporosis ##

outcome_WBC <- outcome_traits_other[194:197,]
outcome_RBC <- outcome_traits_other[158:161,]
outcome_osteo <- extract_outcome_data(snps = exposure_CCR5$SNP, outcome = "UKB-a:87" )

## Sensitivity analysis ##

mr_senstivity_analysis <- function(outcome_trait){
  
  dat <- harmonise_data(exposure_CCR5, outcome_trait)
  
  mr <- mr(dat)
  
  mr_het <- mr_heterogeneity(dat, method_list=c("mr_ivw"))
  
  mr_sct_plot <- mr_scatter_plot(mr, dat)
  
  mr_forest_plot <- mr_forest_plot(mr_singlesnp(dat, all_method = "mr_ivw"))
  
  mr_loo <- mr_leaveoneout(dat)
  
  mr_loo_plot <- mr_leaveoneout_plot(mr_leaveoneout(dat))
  
  
  return(list(mr_het, mr_loo, mr_sct_plot, mr_forest_plot, mr_loo_plot))
  
}

out_WBC <- mr_senstivity_analysis(outcome_WBC)
out_RBC <- mr_senstivity_analysis(outcome_RBC)
out_osteo <- mr_senstivity_analysis(outcome_osteo)

write.table(ldply(out_WBC[1:2], data.frame), file="WBC_het_loo.txt", quote=F, sep="\t")
write.table(ldply(out_RBC[1:2], data.frame), file="RBC_het_loo.txt", quote=F, sep="\t")
write.table(ldply(out_osteo[1:2], data.frame), file="osteo_het_loo.txt", quote=F, sep="\t")




