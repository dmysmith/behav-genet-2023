################################

# Create csv of phenotypes of interest, including residualization
# Diana Smith
# April 2022

rm(list=ls())
################################
# The following R packages need to be loaded

#library(tidyverse)
#library(psych)
library(plyr)
library(dplyr)
#library(PerformanceAnalytics)
#library(pracma)

################################
# This section defines input and output paths for files and functions called. 

# Define the path to the directory which contains the tabulated ABCD data 
inpath <- '/space/abcd-sync/4.0/tabulated/released'
supportpath <- '/space/abcd-sync/4.0/support_files/ABCD_rel4.0_unfiltered'

# Define the full path to the output RDS file 
outpath <- '/space/syn50/1/data/ABCD/d9smith/random_effects/behavioral/data/pheno'

################################
# This code uses helper functions defined in the cmig_tools repo.
# To clone the repo, use:
# git clone https://github.com/cmig-research-group/cmig_tools.git

funcpath <- '/home/d9smith/github/cmig_tools/cmig_tools_utils/r'

source(paste0(funcpath, '/', 'loadtxt.R'))
source(paste0(funcpath, '/', 'makeDEAPdemos.R'))

# The functionmakeDEAPdemos.R requires the path to the directory which 
# contains the tabulated ABCD data defined explicitly here
deappath <- inpath

# Define the file names for the instruments from which we need to pull
# variables. 
# tbx_file <- 'abcd_tbss01.txt'
lmt_file <- 'lmtp201.txt'
# reasoning_file <- 'abcd_ps01.txt'
support_file = 'ABCD_rel4.0_allvars_base_2yr.txt'
phys_file <- 'abcd_ant01.txt'
# Define the full paths to these files 
support_file <- paste0(supportpath, '/', support_file)
lmt_file <- paste0(inpath, '/', lmt_file)
# reasoning_file <- paste0(inpath, '/', reasoning_file)
phys_file <- paste0(inpath, '/', phys_file)

################################
# Load the support file
df <- loadtxt(support_file)

supportvars <- c("src_subject_id","eventname","rel_family_id","interview_age",
"genesis_PC1","genesis_PC2","genesis_PC3","genesis_PC4","genesis_PC5","genesis_PC6","genesis_PC7","genesis_PC8","genesis_PC9","genesis_PC10",
"sex","high.educ","household.income","hisp","race.4level","abcd_site",
"nihtbx_reading_uncorrected","nihtbx_flanker_uncorrected","nihtbx_cardsort_uncorrected","nihtbx_pattern_uncorrected","nihtbx_picture_uncorrected",
"nihtbx_picvocab_uncorrected","nihtbx_list_uncorrected","nihtbx_totalcomp_uncorrected","nihtbx_fluidcomp_uncorrected","nihtbx_cryst_uncorrected",
"pea_ravlt_sd_trial_i_tc","pea_ravlt_sd_trial_ii_tc","pea_ravlt_sd_trial_iii_tc","pea_ravlt_sd_trial_iv_tc","pea_ravlt_sd_trial_v_tc",
"pea_ravlt_sd_trial_sum5trials","pea_ravlt_sd_listb_tc","pea_ravlt_ld_trial_vii_tc","pea_ravlt_learning_slope","pea_wiscv_trs")
fullmat <- df[,supportvars]

################################
# Load the physical variables file  
phys <- loadtxt(phys_file)
# Extract the variables of interest
physvar <- c('src_subject_id', 'eventname', 'anthroheightcalc')
phys <- phys[,physvar]
# Combine with fullmat. 
fullmat <- join(fullmat, phys, by=c('src_subject_id', 'eventname'))

################################
# Load the little man task file
lmt <- loadtxt(lmt_file)
# Extract the variables of interest
lmtvar <- c('src_subject_id', 'eventname', 'lmt_scr_perc_correct')
lmt <- lmt[,lmtvar]
# Combine with the physical health variables. 
fullmat <- join(fullmat, lmt, by=c('src_subject_id', 'eventname'))

################################
# Omit participants with height < 20 or > 80

fullmat <- fullmat[fullmat$anthroheightcalc >= 20 & fullmat$anthroheightcalc <= 80,]

################################
# Create dataframe "baseline" that includes all variables at baseline
baselinevars <- c( "pea_wiscv_trs", "nihtbx_pattern_uncorrected", "nihtbx_flanker_uncorrected", 
        "nihtbx_cardsort_uncorrected", "nihtbx_list_uncorrected", "nihtbx_picture_uncorrected", 
        "nihtbx_picvocab_uncorrected", "nihtbx_reading_uncorrected", "nihtbx_cryst_uncorrected",
        "nihtbx_fluidcomp_uncorrected", "nihtbx_totalcomp_uncorrected", "anthroheightcalc")        

baseline <- fullmat[fullmat$eventname=='baseline_year_1_arm_1',c('src_subject_id','eventname',baselinevars)]


################################
# Create dataframe "longitudinal" that includes baseline and year 2 for all variables with data
y2vars <- c('nihtbx_pattern_uncorrected','nihtbx_flanker_uncorrected','nihtbx_picture_uncorrected',
'nihtbx_picvocab_uncorrected','nihtbx_reading_uncorrected','nihtbx_cryst_uncorrected','anthroheightcalc')

longitudinal <- fullmat[,c('src_subject_id','eventname',y2vars)]

################################
# Save both "baseline" and "longitudinal" as .txt files -- these are unadjusted which we do not use in the manuscript.

if ( ! dir.exists(outpath) ) {
        dir.create(outpath, recursive=TRUE)
}

write.table(baseline, file=paste0(outpath, '/', 'baseline_unadjusted.txt'), sep = "\t", row.names = FALSE)
write.table(longitudinal, file=paste0(outpath, '/', 'longitudinal_unadjusted.txt'), sep = "\t", row.names = FALSE)

## Save variable names to be read by plotting jupyter notebook
write.table(baselinevars, file=paste0(outpath, '/', 'baseline_phenonames.txt'), sep = "\t", row.names = FALSE)
write.table(y2vars, file=paste0(outpath, '/', 'longitudinal_phenonames.txt'), sep = "\t", row.names = FALSE)

################################
# Create several residualized .csv files
################################

# Get just the first 10 genetic PCs and write to a dataframe  
PCs <- c('PC1','PC2','PC3','PC4','PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10')
# Combine with the physical health variables. 
# fullmat <- join(fullmat, pc, by='src_subject_id', match = "all")
fullmat_untouched = fullmat
fullmat = fullmat[,!names(fullmat) %in% c("hisp","race.4level")]

fullmat <- fullmat %>% rename_at(vars(starts_with('genesis_')), ~ PCs)

covariates = c("interview_age",PCs,"sex","high.educ","household.income","abcd_site")

baseline_full = na.omit(fullmat[fullmat$eventname=='baseline_year_1_arm_1',c("src_subject_id","eventname","rel_family_id", covariates, baselinevars)])
longitudinal_full = na.omit(fullmat[,c("src_subject_id","eventname","rel_family_id", covariates, y2vars)])

# create practice effect var for each y2 task
longitudinal_full$prac = 1
dup_ids = longitudinal_full[duplicated(longitudinal_full$src_subject_id),'src_subject_id']
idx = which(longitudinal_full$eventname=='2_year_follow_up_y_arm_1' & longitudinal_full$src_subject_id %in% dup_ids)
longitudinal_full[idx,]$prac = 2

# 1. baseline_full_res_agesexsite
baseline_full_res_agesexsite = baseline_full[,c("src_subject_id", "eventname", baselinevars)]
allModelsList <- lapply(paste(baselinevars, "~ interview_age + sex + abcd_site"), as.formula)
allModelsResults <- lapply(allModelsList, function(x) lm(x, data = baseline_full, na.action = na.exclude))
allModelsResiduals <- lapply(allModelsList, function(x) residuals(lm(x, data = baseline_full)))  
baseline_full_res_agesexsite[,-(1:2)] = allModelsResiduals

# 1b. baseline_full_res_agesex
baseline_full_res_agesex = baseline_full[,c("src_subject_id", "eventname", baselinevars)]
allModelsList <- lapply(paste(baselinevars, "~ interview_age + sex"), as.formula)
allModelsResults <- lapply(allModelsList, function(x) lm(x, data = baseline_full, na.action = na.exclude))
allModelsResiduals <- lapply(allModelsList, function(x) residuals(lm(x, data = baseline_full)))  
baseline_full_res_agesex[,-(1:2)] = allModelsResiduals

# 2. baseline_twins_res_agesexsite
twin_ids = loadtxt(file="/home/d9smith/projects/random_effects/behavioral/twinfiles/twin_IDs.txt")

# check whether both twins have complete data
twin_ids$IID1_complete = twin_ids$IID1 %in% baseline_full$src_subject_id
twin_ids$IID2_complete = twin_ids$IID2 %in% baseline_full$src_subject_id

twin_complete = twin_ids[twin_ids$IID1_complete==T & twin_ids$IID2_complete==T,]
twin_id_list = c(as.character(twin_complete[,1]), as.character(twin_complete[,2]))

tmp = !duplicated(twin_id_list)
IID1_include = tmp[1:dim(twin_complete)[1]]
IID2_include = tmp[(dim(twin_complete)[1]+1):length(tmp)]

twin_unique = twin_complete[IID1_include & IID2_include,]
twin_id_unique = c(as.character(twin_unique[,1]), as.character(twin_unique[,2])) 
twinmat = baseline_full[baseline_full$src_subject_id %in% twin_id_unique,]

baseline_twins_res_agesexsite = twinmat[, c("src_subject_id", "eventname", baselinevars)]
allModelsList <- lapply(paste(baselinevars, "~ interview_age + sex + abcd_site"), as.formula)
allModelsResults <- lapply(allModelsList, function(x) lm(x, data = twinmat, na.action = na.exclude))
allModelsResiduals <- lapply(allModelsList, function(x) residuals(lm(x, data = twinmat)))  
baseline_twins_res_agesexsite[,-(1:2)] = allModelsResiduals

# 2b. baseline_twins_res_agesex 
baseline_twins_res_agesex = twinmat[, c("src_subject_id", "eventname", baselinevars)]
allModelsList <- lapply(paste(baselinevars, "~ interview_age + sex"), as.formula)
allModelsResults <- lapply(allModelsList, function(x) lm(x, data = twinmat, na.action = na.exclude))
allModelsResiduals <- lapply(allModelsList, function(x) residuals(lm(x, data = twinmat)))  
baseline_twins_res_agesex[,-(1:2)] = allModelsResiduals

# save file with completeness info
write.table(twin_unique, file="/home/d9smith/projects/random_effects/behavioral/twinfiles/twin_IDs_complete.txt", sep = "\t", row.names = FALSE)
write.table(twin_unique, file="/home/d9smith/projects/random_effects/behavioral/twinfiles/twin_IDs_complete.txt", sep = "\t", row.names = FALSE)

# 3. baseline_full_res_agesexsitepcs
baseline_full_res_agesexsitepcs = baseline_full[,c("src_subject_id", "eventname", baselinevars)]
allModelsList <- lapply(paste(baselinevars, "~ interview_age + sex + abcd_site + PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10"), as.formula)
allModelsResults <- lapply(allModelsList, function(x) lm(x, data = baseline_full, na.action = na.exclude))
allModelsResiduals <- lapply(allModelsList, function(x) residuals(lm(x, data = baseline_full)))  
baseline_full_res_agesexsitepcs[,-(1:2)] = allModelsResiduals

# 4. baseline_full_res_agesexsiteeducinc
baseline_full_res_agesexsiteeducinc = baseline_full[,c("src_subject_id", "eventname", baselinevars)]
allModelsList <- lapply(paste(baselinevars, "~ interview_age + sex + abcd_site + high.educ + household.income"), as.formula)
allModelsResults <- lapply(allModelsList, function(x) lm(x, data = baseline_full, na.action = na.exclude))
allModelsResiduals <- lapply(allModelsList, function(x) residuals(lm(x, data = baseline_full)))  
baseline_full_res_agesexsiteeducinc[,-(1:2)] = allModelsResiduals

# 5. baseline_full_res_agesexeducincpcs
baseline_full_res_agesexsiteeducincpcs = baseline_full[,c("src_subject_id", "eventname", baselinevars)]
allModelsList <- lapply(paste(baselinevars, "~ interview_age + sex + abcd_site + high.educ + household.income + PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10"), as.formula)
allModelsResults <- lapply(allModelsList, function(x) lm(x, data = baseline_full, na.action = na.exclude))
allModelsResiduals <- lapply(allModelsList, function(x) residuals(lm(x, data = baseline_full)))  
baseline_full_res_agesexsiteeducincpcs[,-(1:2)] = allModelsResiduals

# 6. twins_res_agesexsiteeducincpcs
baseline_twins_res_agesexsiteeducincpcs = twinmat[,c("src_subject_id", "eventname", baselinevars)]
allModelsList <- lapply(paste(baselinevars, "~ interview_age + sex + abcd_site + high.educ + household.income + PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10"), as.formula)
allModelsResults <- lapply(allModelsList, function(x) lm(x, data = twinmat, na.action = na.exclude))
allModelsResiduals <- lapply(allModelsList, function(x) residuals(lm(x, data = twinmat)))  
baseline_twins_res_agesexsiteeducincpcs[,-(1:2)] = allModelsResiduals 

# 7. longitudinal_full_res_agesexsiteprac
longitudinal_full_res_agesexsiteprac = longitudinal_full[,c("src_subject_id", "eventname", y2vars)]
allModelsList <- lapply(paste(y2vars, "~ interview_age + sex + abcd_site + prac"), as.formula)
allModelsResults <- lapply(allModelsList, function(x) lm(x, data = longitudinal_full, na.action = na.exclude))
allModelsResiduals <- lapply(allModelsList, function(x) residuals(lm(x, data = longitudinal_full)))  
longitudinal_full_res_agesexsiteprac[,-(1:2)] = allModelsResiduals

# 7b. longitudinal_full_res_agesexprac
longitudinal_full_res_agesexprac = longitudinal_full[,c("src_subject_id", "eventname", y2vars)]
allModelsList <- lapply(paste(y2vars, "~ interview_age + sex + prac"), as.formula)
allModelsResults <- lapply(allModelsList, function(x) lm(x, data = longitudinal_full, na.action = na.exclude))
allModelsResiduals <- lapply(allModelsList, function(x) residuals(lm(x, data = longitudinal_full)))  
longitudinal_full_res_agesexprac[,-(1:2)] = allModelsResiduals 

# 8. longitudinal_full_res_agesexsitepraceducincpcs
longitudinal_full_res_agesexsitepraceducincpcs = longitudinal_full[,c("src_subject_id", "eventname", y2vars)]
allModelsList <- lapply(paste(y2vars, "~ interview_age + sex + abcd_site + prac + high.educ + household.income + PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10"), as.formula)
allModelsResults <- lapply(allModelsList, function(x) lm(x, data = longitudinal_full, na.action = na.exclude))
allModelsResiduals <- lapply(allModelsList, function(x) residuals(lm(x, data = longitudinal_full)))  
longitudinal_full_res_agesexsitepraceducincpcs[,-(1:2)] = allModelsResiduals 

# 9. longitudinal_notwins_res_agesexsiteprac
all_twin_ids = c(as.character(twin_ids$IID1),as.character(twin_ids$IID2))
longitudinal_notwins = longitudinal_full[-which(longitudinal_full$src_subject_id %in% all_twin_ids),]

longitudinal_notwins_res_agesexsiteprac = longitudinal_notwins[,c("src_subject_id", "eventname", y2vars)]
allModelsList <- lapply(paste(y2vars, "~ interview_age + sex + abcd_site + prac"), as.formula)
allModelsResults <- lapply(allModelsList, function(x) lm(x, data = longitudinal_notwins, na.action = na.exclude))
allModelsResiduals <- lapply(allModelsList, function(x) residuals(lm(x, data = longitudinal_notwins)))  
longitudinal_notwins_res_agesexsiteprac[,-(1:2)] = allModelsResiduals

# 10. longitudinal_full_res_agesexsite
longitudinal_full_res_agesexsite = longitudinal_full[,c("src_subject_id", "eventname", y2vars)]
allModelsList <- lapply(paste(y2vars, "~ interview_age + sex + abcd_site"), as.formula)
allModelsResults <- lapply(allModelsList, function(x) lm(x, data = longitudinal_full, na.action = na.exclude))
allModelsResiduals <- lapply(allModelsList, function(x) residuals(lm(x, data = longitudinal_full)))  
longitudinal_full_res_agesexsite[,-(1:2)] = allModelsResiduals 

# 11. longitudinal_full_res_agesexsiteeducincpcs
longitudinal_full_res_agesexsiteeducincpcs = longitudinal_full[,c("src_subject_id", "eventname", y2vars)]
allModelsList <- lapply(paste(y2vars, "~ interview_age + sex + abcd_site + high.educ + household.income + PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10"), as.formula)
allModelsResults <- lapply(allModelsList, function(x) lm(x, data = longitudinal_full, na.action = na.exclude))
allModelsResiduals <- lapply(allModelsList, function(x) residuals(lm(x, data = longitudinal_full)))  
longitudinal_full_res_agesexsiteeducincpcs[,-(1:2)] = allModelsResiduals 

# 12. longitudinal_notwins_res_agesexsite
all_twin_ids = c(as.character(twin_ids$IID1),as.character(twin_ids$IID2))
longitudinal_notwins = longitudinal_full[-which(longitudinal_full$src_subject_id %in% all_twin_ids),]

longitudinal_notwins_res_agesexsite = longitudinal_notwins[,c("src_subject_id", "eventname", y2vars)]
allModelsList <- lapply(paste(y2vars, "~ interview_age + sex + abcd_site"), as.formula)
allModelsResults <- lapply(allModelsList, function(x) lm(x, data = longitudinal_notwins, na.action = na.exclude))
allModelsResiduals <- lapply(allModelsList, function(x) residuals(lm(x, data = longitudinal_notwins)))  
longitudinal_notwins_res_agesexsite[,-(1:2)] = allModelsResiduals

# 12b. longitudinal_notwins_res_agesexprac 
longitudinal_notwins_res_agesexprac = longitudinal_notwins[,c("src_subject_id", "eventname", y2vars)]
allModelsList <- lapply(paste(y2vars, "~ interview_age + sex + prac"), as.formula)
allModelsResults <- lapply(allModelsList, function(x) lm(x, data = longitudinal_notwins, na.action = na.exclude))
allModelsResiduals <- lapply(allModelsList, function(x) residuals(lm(x, data = longitudinal_notwins)))  
longitudinal_notwins_res_agesexprac[,-(1:2)] = allModelsResiduals

# 13. y2_full_res_agesexsite - we did not use this in the paper
y2_full = longitudinal_full[longitudinal_full$eventname=="2_year_follow_up_y_arm_1",]

y2_full_res_agesexsite = y2_full[,c("src_subject_id", "eventname", y2vars)]
allModelsList <- lapply(paste(y2vars, "~ interview_age + sex + abcd_site"), as.formula)
allModelsResults <- lapply(allModelsList, function(x) lm(x, data = y2_full, na.action = na.exclude))
allModelsResiduals <- lapply(allModelsList, function(x) residuals(lm(x, data = y2_full)))  
y2_full_res_agesexsite[,-(1:2)] = allModelsResiduals

# 13b. y2_full_res_agesex - not used in the paper
y2_full_res_agesex = y2_full[,c("src_subject_id", "eventname", y2vars)]
allModelsList <- lapply(paste(y2vars, "~ interview_age + sex"), as.formula)
allModelsResults <- lapply(allModelsList, function(x) lm(x, data = y2_full, na.action = na.exclude))
allModelsResiduals <- lapply(allModelsList, function(x) residuals(lm(x, data = y2_full)))  
y2_full_res_agesex[,-(1:2)] = allModelsResiduals

# 14. y2_full_res_agesexsiteeducincpcs - not used in the paper
y2_full_res_agesexsiteeducincpcs = y2_full[,c("src_subject_id", "eventname", y2vars)]
allModelsList <- lapply(paste(y2vars, "~ interview_age + sex + abcd_site + high.educ + household.income + PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10"), as.formula)
allModelsResults <- lapply(allModelsList, function(x) lm(x, data = y2_full, na.action = na.exclude))
allModelsResiduals <- lapply(allModelsList, function(x) residuals(lm(x, data = y2_full)))  
y2_full_res_agesexsiteeducincpcs[,-(1:2)] = allModelsResiduals

# save phenofiles  
write.table(baseline_full_res_agesexsite, file=paste0(outpath, '/', 'baseline_full_res_agesexsite.txt'), sep = "\t", row.names = FALSE)
write.table(baseline_full_res_agesex, file=paste0(outpath, '/', 'baseline_full_res_agesex.txt'), sep = "\t", row.names = FALSE)
write.table(baseline_twins_res_agesexsite, file=paste0(outpath, '/', 'baseline_twins_res_agesexsite.txt'), sep = "\t", row.names = FALSE)
write.table(baseline_twins_res_agesex, file=paste0(outpath, '/', 'baseline_twins_res_agesex.txt'), sep = "\t", row.names = FALSE)
write.table(baseline_full_res_agesexsitepcs, file=paste0(outpath, '/', 'baseline_full_res_agesexsitepcs.txt'), sep = "\t", row.names = FALSE)
write.table(baseline_full_res_agesexsiteeducinc, file=paste0(outpath, '/', 'baseline_full_res_agesexsiteeducinc.txt'), sep = "\t", row.names = FALSE)
write.table(baseline_full_res_agesexsiteeducincpcs, file=paste0(outpath, '/', 'baseline_full_res_agesexsiteeducincpcs.txt'), sep = "\t", row.names = FALSE)
write.table(baseline_twins_res_agesexsiteeducincpcs, file=paste0(outpath, '/', 'baseline_twins_res_agesexsiteeducincpcs.txt'), sep = "\t", row.names = FALSE)

write.table(longitudinal_full_res_agesexsiteprac, file=paste0(outpath, '/', 'longitudinal_full_res_agesexsiteprac.txt'), sep = "\t", row.names = FALSE)
write.table(longitudinal_full_res_agesexprac, file=paste0(outpath, '/', 'longitudinal_full_res_agesexprac.txt'), sep = "\t", row.names = FALSE)
write.table(longitudinal_full_res_agesexsitepraceducincpcs, file=paste0(outpath, '/', 'longitudinal_full_res_agesexsitepraceducincpcs.txt'), sep = "\t", row.names = FALSE)
write.table(longitudinal_notwins_res_agesexsiteprac, file=paste0(outpath, '/', 'longitudinal_notwins_res_agesexsiteprac.txt'), sep = "\t", row.names = FALSE)

write.table(longitudinal_full_res_agesexsite, file=paste0(outpath, '/', 'longitudinal_full_res_agesexsite.txt'), sep = "\t", row.names = FALSE)
write.table(longitudinal_full_res_agesexsiteeducincpcs, file=paste0(outpath, '/', 'longitudinal_full_res_agesexsiteeducincpcs.txt'), sep = "\t", row.names = FALSE)
write.table(longitudinal_notwins_res_agesexsite, file=paste0(outpath, '/', 'longitudinal_notwins_res_agesexsite.txt'), sep = "\t", row.names = FALSE)
write.table(longitudinal_notwins_res_agesexprac, file=paste0(outpath, '/', 'longitudinal_notwins_res_agesexprac.txt'), sep = "\t", row.names = FALSE)

write.table(y2_full_res_agesexsite, file=paste0(outpath, '/', 'y2_full_res_agesexsite.txt'), sep = "\t", row.names = FALSE)
write.table(y2_full_res_agesex, file=paste0(outpath, '/', 'y2_full_res_agesex.txt'), sep = "\t", row.names = FALSE)
write.table(y2_full_res_agesexsiteeducincpcs, file=paste0(outpath, '/', 'y2_full_res_agesexsiteeducincpcs.txt'), sep = "\t", row.names = FALSE)