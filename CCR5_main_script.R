setwd("/Users/hw15842/Documents/PhD work/CCR5_paper")

library (devtools)
library(plyr)

install_github("MRCIEU/TwoSampleMR")

library(TwoSampleMR)

ao <- available_outcomes()

CCR5_eQTL_gen <- read.table("CCR5.txt", header=TRUE)

outcome_ids <- read.csv("trait_list_for_Hannah.csv", header=F)

run_mr_CCR5 <- function(outcome_ID){
  
  exposure_1 <- CCR5_eQTL_gen
  names(exposure_1) <- c("SNP", "Phenotype", "id.exposure", "CHR",  "BP", "effect_allele.exposure","other_allele.exposure", 
                         "eaf.exposure", "beta.exposure", "se.exposure", "pval.exposure")
  exposure_1$exposure <- "CCR5_ENSG00000160791"
  
  outcome_1 <- extract_outcome_data(snps = exposure_1$SNP, outcome = outcome_ID)
  dat1 <- harmonise_data(exposure_1, outcome_1)
  res1 <- mr(dat1, method_list = "mr_ivw")
  return(res1)
}

outcome_ids_rowsremoved <- outcome_ids[-c(548,650,709, 710, 758, 759, 760, 761),] # removed rows that wont run as didnt have the CCR5 snps
results_all_rowsremoved <- lapply(outcome_ids_rowsremoved$V1, run_mr_CCR5)

first250 <- head(outcome_ids_rowsremoved, n=250)
results_first250 <- lapply(first250$V1, run_mr_CCR5)
results_first250 <- ldply(results_first250, data.frame)

out251to400 <- outcome_ids_rowsremoved[251:400,]
results_251to400 <- lapply(out251to400$V1, run_mr_CCR5)
results_251to400 <- ldply(results_251to400, data.frame)

out401to600 <- outcome_ids_rowsremoved[401:600,]
results_401to600 <- lapply(out401to600$V1, run_mr_CCR5)
results_401to600 <- ldply(results_401to600, data.frame)

out601to674 <- outcome_ids_rowsremoved[601:674,]
results_601to674 <- lapply(out601to674$V1, run_mr_CCR5)
results_601to674 <- ldply(results_601to674, data.frame)

out675to754 <- outcome_ids_rowsremoved[675:754,]
results_675to754 <- lapply(out675to754$V1, run_mr_CCR5)
results_675to754 <- ldply(results_675to754, data.frame)

results_all <- rbind(results_first250, results_251to400, results_401to600, results_601to674, results_675to754)


results_all_rowsreomved_dataframe <- ldply(results_all_rowsremoved, data.frame)

write.table(results_all, file="results_all_tabsep.txt", quote=F, sep="\t")

results_all_read_in <- read.table("results_all_tabsep.txt", header=T, sep="\t")



## Sort by pval ##

results_sorted_pval <- results_all[order(results_all$pval),]

write.table(results_sorted_pval, file = "results_sorted_pval.txt", quote=F, sep="\t")

0.05/674




#### Now need to look at just the cognitive traits...###

head(outcome_ids)

cognitive_outcomes <- subset(outcome_ids, grepl("reaction time|memory|education|reasoning", V2)) ## these were from other studies that have looked at cognitive in biobank

results_cognitive <- lapply(cognitive_outcomes$V1, run_mr_CCR5)
results_cognitive <- ldply(results_cognitive, data.frame)
write.table(results_cognitive, file="results_cognitive.txt", quote=F, sep="\t")

## Only education and memory have results and neither of those is significant... ##

## All Diseases ##

disease_outcomes <- subset(outcome_ids, grepl("disease", V2))
disease_outcomes <- disease_outcomes[-c(28),]
results_disease <- lapply(disease_outcomes$V1, run_mr_CCR5)
results_disease <- ldply(results_disease, data.frame)
write.table(results_disease, file="results_disease.txt", quote=F, sep="\t")

## Ceoliac disease only one that looks a little bit significant and that has been shown in the literature that CCR5delta32 increases risk of ceolic 

0.05/22

###########################################
## To read the tables back in if need be ##
###########################################

results_all <- read.table("results_all_tabsep.txt", header=T, sep="\t")
results_sorted_pval <- read.table("results_sorted_pval.txt", header=T, sep="\t")
results_cognitive <- read.table("results_cognitive.txt", header=T, sep="\t")
results_disease <- read.table("results_disease.txt", header=T, sep="\t")


### Sensitivity analysis for Osteoporosis (nearly significant in results_all) ##


exposure_1 <- CCR5_eQTL_gen
names(exposure_1) <- c("SNP", "Phenotype", "id.exposure", "CHR",  "BP", "effect_allele.exposure","other_allele.exposure", 
                       "eaf.exposure", "beta.exposure", "se.exposure", "pval.exposure")
exposure_1$exposure <- "CCR5_ENSG00000160791" 

outcome_1 <- extract_outcome_data(snps = exposure_1$SNP, outcome = "UKB-a:87" )

dat_osteo <- harmonise_data(exposure_1, outcome_1)

osteo_mr_res <- mr(dat_osteo)

osteo_mr_het <- mr_heterogeneity(dat_osteo)

osteo_mr_sct_plot <- mr_scatter_plot(osteo_mr_res, dat_osteo)

osteo_mr_forest_plot <- mr_forest_plot(mr_singlesnp(dat_osteo, all_method="mr_ivw"))

osteo_mr_loo <- mr_leaveoneout(dat_osteo)
osteo_mr_loo_plot <-mr_leaveoneout_plot(mr_leaveoneout(dat_osteo))

osteo_mr_plei <- mr_pleiotropy_test(dat_osteo)



### Sensitivity analysis for Celiac disease ###

outcome_2 <- extract_outcome_data(snps = exposure_1$SNP, outcome = "UKB-a:101" )

dat_celiac <- harmonise_data(exposure_1, outcome_2)

celiac_mr_res <- mr(dat_celiac)

celiac_mr_het <- mr_heterogeneity(dat_celiac)

celiac_mr_sct_plot <- mr_scatter_plot(celiac_mr_res, dat_celiac)

celiac_mr_forest_plot <- mr_forest_plot(mr_singlesnp(dat_celiac, all_method="mr_ivw"))

celiac_mr_loo <- mr_leaveoneout(dat_celiac)
celiac_mr_loo_plot <- mr_leaveoneout_plot(mr_leaveoneout(dat_celiac))

celiac_mr_plei <-mr_pleiotropy_test(dat_celiac)












## Perfrom Radial MR? ###

install.packages("devtools")
library(devtools)
install_github("WSpiller/RadialMR")
