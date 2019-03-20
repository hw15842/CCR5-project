setwd("/Users/hw15842/Documents/PhD work/CCR5_paper")

library (devtools)

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

outcome_ids_rowsremoved <- outcome_ids[-c(548,650,709, 710, 758, 759, 760, 761),]
results_all_rowsremoved <- lapply(outcome_ids_rowsremoved$V1, run_mr_CCR5)


results_all_rowsreomved_dataframe <- ldply(results_all_rowsremoved, data.frame)

write.table(results_all_rowsreomved_dataframe, file="results_all_rowsreomved_dataframe.txt", quote=F)

