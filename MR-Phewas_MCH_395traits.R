##############################################################
### MR-Phewas of mean corpuscular hemoglobin on 395 traits ###
##############################################################




setwd("/newhome/hw15842/PhD/CCR5/MR-Phewas_MCH")

install.packages("devtools", lib = "/newhome/hw15842/R/x86_64-pc-linux-gnu-library/3.5", repos='http://cran.us.r-project.org')
install.packages("plyr", lib = "/newhome/hw15842/R/x86_64-pc-linux-gnu-library/3.5", repos='http://cran.us.r-project.org')
install.packages("data.table", lib = "/newhome/hw15842/R/x86_64-pc-linux-gnu-library/3.5", repos='http://cran.us.r-project.org')

library (devtools)
library(plyr)
library(data.table)

install_github("MRCIEU/TwoSampleMR")

library(TwoSampleMR)

ao <- available_outcomes()


## exposure data ##

MCH_all_data <- fread("blood_MEAN_CORPUSCULAR_HEMOGLOBIN.sumstats")

names(MCH_all_data) <- c("SNP", "CHR", "position", "effect_allele", "other_allele", "reference_allele", "eaf", "beta", "se", "pval", "N", "info")

matching <- match(MCH_all_data$effect_allele, MCH_all_data$other_allele)

sum(is.na(matching))

MCH_all_data$MAF <- ifelse (MCH_all_data$eaf<=0.5, as.numeric(MCH_all_data$eaf),
                    ifelse(MCH_all_data$eaf>=0.5, as.numeric(1-(MCH_all_data$eaf)), NA))



MCH_all_data <- MCH_all_data[!(MCH_all_data$MAF<=0.001 & MCH_all_data$info<=0.8),]

exposure_MCH <- format_data(MCH_all_data, type="exposure", snp_col="SNP", beta_col="beta", se_col="se", eaf_col="eaf", effect_allele_col="effect_allele", other_allele_col="other_allele", pval_col="pval", pos_col="position", chr_col="CHR")


head(exposure_MCH)
nrow(exposure_MCH)

## CLUMP

exposure_MCH <- clump_data(exposure_MCH)

## MRbase outcome traits ##

outcome_traits_mrbase <- read.csv("list_of_395_traits.csv", header=T)

#remove NAs 
outcome_traits_mrbase <- outcome_traits_mrbase[!is.na(outcome_traits_mrbase$MR.Base.id),]




run_mr_MCH <- function(outcome_ID){
  
  outcome_1 <- extract_outcome_data(snps = exposure_MCH$SNP, outcome = outcome_ID)
  dat1 <- harmonise_data(exposure_MCH, outcome_1)
  res1 <- mr(dat1, method_list = "mr_ivw")
  return(res1)
}

results_MRbase_MCH <- lapply(outcome_traits_mrbase$MR.Base.id, run_mr_MCH)
results_MRbase_MCH_table <- ldply(results_MRbase_MCH, data.frame)



write.table(results_MRbase__MCH_table, file = "results_MRbase_MCH.txt", quote=F, sep="\t")



