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



#MCH_all_data <- MCH_all_data[!(MCH_all_data$MAF<=0.001 & MCH_all_data$info<=0.8 & MCH_all_data$pval>=5e-8),]

nrow(MCH_all_data)

MCH_all_data <- MCH_all_data[which(MCH_all_data$MAF>=0.001 & MCH_all_data$info>=0.8 & MCH_all_data$pval<=5e-8),]

nrow(MCH_all_data)

exposure_MCH <- format_data(MCH_all_data, type="exposure", snp_col="SNP", beta_col="beta", se_col="se", eaf_col="eaf", effect_allele_col="effect_allele", other_allele_col="other_allele", pval_col="pval", pos_col="position", chr_col="CHR")

nrow(exposure_MCH)


## CLUMP


data_split <- split(exposure_MCH, (as.numeric(rownames(exposure_MCH))-1) %/% 1000)


data_by_chunks <- function(section){
  clumped_data_out <- suppress_messages((clump_data(data_split[[section]])), pattern="rs")
  savefile1=paste("clumped_data_section",section, sep="_")
  #saveRDS(clumped_data_out, file = savefile1)
  assign(savefile1, clumped_data_out,envir = globalenv())
  return(clumped_data_out)
}


section_numbers <- c(1:(length(data_split)))

try(exposure_MCH_split <- lapply(section_numbers, data_by_chunks))

#try(exposure_MCH_split <- ldply(exposure_MCH_split, data.frame))

ls()

results_all <- rbindlist(mget(ls( pattern = "clumped_data_section")))

head(results_all)
nrow(results_all)

ls()

##### "try" function should allow the rest of the script to run, now need something in here that will say, "where did it get to first time, pick up from there and do it agin until finished"

a<-length(ls(pattern="clumped_data_section"))


#section_numbers <- c(length(ls(pattern="clumped_data_section"))):(length(data_split)))
#try(exposure_MCH_split <- lapply(section_numbers, data_by_chunks))



i=(length(ls(pattern="clumped_data_section")))
   while( i < (length(data_split))){
     
     if ((length(ls(pattern="clumped_data_section")))<(length(data_split))) { 
       section_numbers <- c(length(ls(pattern="clumped_data_section"))):(length(data_split))
		try(exposure_MCH_split <- lapply(section_numbers, data_by_chunks))
     }else{
       break
     }
   }

ls()

results_all <- rbindlist(mget(ls( pattern = "clumped_data_section")))

head(results_all)
nrow(results_all)

ls()

write.table(results_all, file="clumping_results_all.txt", quote=F, sep="\t")


#section_numbers_split <- split(section_numbers, ceiling(seq_along(section_numbers)/500000))##
#

#section1 <- section_numbers_split[[1]]
#section2 <- section_numbers_split[[2]]
#section3 <- section_numbers_split[[3]]
#section4 <- section_numbers_split[[4]]
#section5 <- section_numbers_split[[5]]
#section6 <- section_numbers_split[[6]]
#section7 <- section_numbers_split[[7]]
#section8 <- section_numbers_split[[8]]
#section9 <- section_numbers_split[[9]]
#section10 <- section_numbers_split[[10]]
#section11 <- section_numbers_split[[11]]
#section12 <- section_numbers_split[[12]]
#section13 <- section_numbers_split[[13]]
#section14 <- section_numbers_split[[14]]


#results1 <- lapply(section1, data_by_chunks)
#results1 <-  ldply(results1, data.frame)
#head(results1)
#nrow(results1)#

#results2 <- lapply(section2, data_by_chunks)
#results2 <-  ldply(results2, data.frame)
#head(results2)
#nrow(results2)#

#results3 <- lapply(section3, data_by_chunks)
#results3 <-  ldply(results3, data.frame)
#results4 <- lapply(section4, data_by_chunks)
#results4 <-  ldply(results4, data.frame)
#results5 <- lapply(section5, data_by_chunks)
#results5 <-  ldply(results5, data.frame)
#results6 <- lapply(section6, data_by_chunks)
#results6 <-  ldply(results6, data.frame)
#results7 <- lapply(section7, data_by_chunks)
#results7 <-  ldply(results7, data.frame)
#results8 <- lapply(section8, data_by_chunks)
#results8 <-  ldply(results8, data.frame)#
#

#head(results1)


#write.table(exposure_MCH_split_first_half, file="exposure_MCH_split_first_half.txt", quote=F, sep="\t")
#write.table(exposure_MCH_split_second_half, file="exposure_MCH_split_second_half.txt", quote=F, sep="\t")





## MRbase outcome traits ##

outcome_traits_mrbase <- read.csv("list_of_395_traits.csv", header=T)

#remove NAs 
outcome_traits_mrbase <- outcome_traits_mrbase[!is.na(outcome_traits_mrbase$MR.Base.id),]




run_mr_MCH <- function(outcome_ID){
  
  outcome_1 <- extract_outcome_data(snps = results_all$SNP, outcome = outcome_ID)
  dat1 <- harmonise_data(results_all, outcome_1)
  res1 <- mr(dat1, method_list = "mr_ivw")
  return(res1)
}

results_MRbase_MCH <- lapply(outcome_traits_mrbase$MR.Base.id, run_mr_MCH)
results_MRbase_MCH_table <- ldply(results_MRbase_MCH, data.frame)



write.table(results_MRbase__MCH_table, file = "results_MRbase_MCH.txt", quote=F, sep="\t")



