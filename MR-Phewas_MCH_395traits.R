##############################################################
### MR-Phewas of mean corpuscular hemoglobin on 395 traits ###
##############################################################




setwd("/newhome/hw15842/PhD/CCR5/MR-Phewas_MCH")

install.packages("devtools", lib = "/newhome/hw15842/R/x86_64-pc-linux-gnu-library/3.5", repos='http://cran.us.r-project.org')
install.packages("plyr", lib = "/newhome/hw15842/R/x86_64-pc-linux-gnu-library/3.5", repos='http://cran.us.r-project.org')
install.packages("data.table", lib = "/newhome/hw15842/R/x86_64-pc-linux-gnu-library/3.5", repos='http://cran.us.r-project.org')
install.packages("pkgcond", lib = "/newhome/hw15842/R/x86_64-pc-linux-gnu-library/3.5", repos='http://cran.us.r-project.org')



library (devtools)
library(plyr)
library(data.table)
library(pkgcond)

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



MCH_all_data <- MCH_all_data[!(MCH_all_data$MAF<=0.001 & MCH_all_data$info<=0.8 & MCH_all_data$pval>=5e-8),]

exposure_MCH <- format_data(MCH_all_data, type="exposure", snp_col="SNP", beta_col="beta", se_col="se", eaf_col="eaf", effect_allele_col="effect_allele", other_allele_col="other_allele", pval_col="pval", pos_col="position", chr_col="CHR")

nrow(exposure_MCH)


## CLUMP


data_split <- split(exposure_MCH, (as.numeric(rownames(exposure_MCH))-1) %/% 1000)

data_by_chunks <- function(section){
  clumped_data_out <- suppress_messages((clump_data(data_split[[section]])), pattern="rs")
  return(clumped_data_out)
}

section_numbers <- c(1:(length(data_split)))

exposure_MCH_split <- lapply(section_numbers, data_by_chunks)

exposure_MCH_split <- ldply(exposure_MCH_split, data.frame)







## MRbase outcome traits ##

outcome_traits_mrbase <- read.csv("list_of_395_traits.csv", header=T)

#remove NAs 
outcome_traits_mrbase <- outcome_traits_mrbase[!is.na(outcome_traits_mrbase$MR.Base.id),]




run_mr_MCH <- function(outcome_ID){
  
  outcome_1 <- extract_outcome_data(snps = exposure_MCH_split$SNP, outcome = outcome_ID)
  dat1 <- harmonise_data(exposure_MCH_split, outcome_1)
  res1 <- mr(dat1, method_list = "mr_ivw")
  return(res1)
}

results_MRbase_MCH <- lapply(outcome_traits_mrbase$MR.Base.id, run_mr_MCH)
results_MRbase_MCH_table <- ldply(results_MRbase_MCH, data.frame)



write.table(results_MRbase__MCH_table, file = "results_MRbase_MCH.txt", quote=F, sep="\t")



