################################

# Create csv of cognition variables to use for FEMA
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
inpath <- '/space/amdale/1/tmp/ABCD_cache/abcd-sync/4.0/tabulated/released'
supportpath <- '/space/amdale/1/tmp/ABCD_cache/abcd-sync/4.0/support_files/ABCD_rel4.0_unfiltered'

# Define the full path to the output RDS file 
outpath <- '/space/syn50/1/data/ABCD/d9smith/random_effects/behavioral/data/pheno'

# Define the path to tge cmig_tools_utils/r directory
funcpath <- '/home/d9smith/github/cmig_tools_internal/cmig_tools_utils/r'
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

# R needs to parse two functions from cmig_utils/r. One is load.txt
# to load the .txt files and remove unnecessary 2nd row, the other is 
# makeDEAPdemos.R which collates SES variables into the format used by
# DESP. 
source(paste0(funcpath, '/', 'loadtxt.R'))
source(paste0(funcpath, '/', 'makeDEAPdemos.R'))

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
# Create dataframe "baseline" that includes all variables at baseline
baselinevars <- names(fullmat)[-(1:20)]
baseline <- fullmat[fullmat$eventname=='baseline_year_1_arm_1',baselinevars]

################################
# Create dataframe "longitudinal" that includes baseline and year 2 for all variables with data
y2vars <- c('src_subject_id','eventname', 'nihtbx_picvocab_uncorrected','nihtbx_flanker_uncorrected','nihtbx_pattern_uncorrected',
'nihtbx_picture_uncorrected','nihtbx_reading_uncorrected','nihtbx_cryst_uncorrected','lmt_scr_perc_correct','anthroheightcalc')

longitudinal <- fullmat[,y2vars]

################################
# Save both "baseline" and "longitudinal" as RDS files 

if ( ! dir.exists(outpath) ) {
        dir.create(outpath, recursive=TRUE)
}

write.table(baseline, file=paste0(outpath, '/', 'baseline_unadjusted.txt'), sep = "\t", row.names = FALSE)
write.table(longitudinal, file=paste0(outpath, '/', 'longitudinal_unadjusted.txt'), sep = "\t", row.names = FALSE)

################################
# Create several residualized .csv files
################################
# Create the SES variables as coded by DEAP
# deap <- makeDEAPdemos(deappath)
# deap <- deap[ , c("src_subject_id", "eventname", "sex", "interview_age", "high.educ", "household.income")]
# Combine with the previously extracted variables
# fullmat <- join(fullmat, deap, by=c('src_subject_id', 'eventname'))

# genetic PCs subject data
# pcfile <- '/space/gwas-syn2/1/data/GWAS/ABCD/genotype_proc/imputation/pop_struct_smokescreen/ABCD_20220428.updated.nodups.curated_pcair.tsv'
# pc_mat <- read.delim(pcfile)
# Get just the first 10 PCs and write to a dataframe  
# pc <- data.frame(pc_mat[,c('X','C1','C2','C3','C4','C5', 'C6', 'C7', 'C8', 'C9', 'C10')])
PCs <- c('PC1','PC2','PC3','PC4','PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10')
# Combine with the physical health variables. 
# fullmat <- join(fullmat, pc, by='src_subject_id', match = "all")
fullmat_untouched = fullmat
fullmat = fullmat[,!names(fullmat) %in% c("hisp","race.4level")]

fullmat <- fullmat %>% rename_at(vars(starts_with('genesis_')), ~ PCs)

covariates = c("interview_age",PCs,"sex","high.educ","household.income","abcd_site")

baseline_full = na.omit(fullmat[fullmat$eventname=='baseline_year_1_arm_1',c("src_subject_id","eventname","rel_family_id", covariates, baselinevars)])
longitudinal_full = na.omit(fullmat[,c("src_subject_id","eventname","rel_family_id", covariates, y2vars)])

# 1. baseline_full_res_agesexsite
baseline_full_res_agesexsite = baseline_full[,c("src_subject_id", "eventname", baselinevars)]
allModelsList <- lapply(paste(baselinevars, "~ interview_age + sex + abcd_site"), as.formula)
allModelsResults <- lapply(allModelsList, function(x) lm(x, data = baseline_full, na.action = na.exclude))
allModelsResiduals <- lapply(allModelsList, function(x) residuals(lm(x, data = baseline_full)))  
baseline_full_res_agesexsite[,-(1:2)] = allModelsResiduals

# 2. baseline_twins_res_agesexsite
twin_ids = loadtxt(file="/home/d9smith/projects/random_effects/behavioral/twinfiles/twin_IDs.txt")

# check whether both twins have complete data
twin_ids$IID1_complete = twin_ids$IID1 %in% baseline_full$src_subject_id
twin_ids$IID2_complete = twin_ids$IID2 %in% baseline_full$src_subject_id
# save file with completeness info
write.table(twin_ids, file="/home/d9smith/projects/random_effects/behavioral/twinfiles/twin_IDs_complete.txt")

twin_complete = twin_ids[twin_ids$IID1_complete==T & twin_ids$IID2_complete==T,]

twin_id_list = c(as.character(twin_complete[,1]), as.character(twin_complete[,2]))
twinmat = baseline_full[baseline_full$src_subject_id %in% twin_id_list,]
baseline_twins_res_agesexsite = twinmat[, c("src_subject_id", "eventname", baselinevars)]

allModelsList <- lapply(paste(baselinevars, "~ interview_age + sex + abcd_site"), as.formula)
allModelsResults <- lapply(allModelsList, function(x) lm(x, data = twinmat, na.action = na.exclude))
allModelsResiduals <- lapply(allModelsList, function(x) residuals(lm(x, data = twinmat)))  
baseline_twins_res_agesexsite[,-(1:2)] = allModelsResiduals

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

# save phenofiles - all baseline for now - DS 2022-08-16
write.table(baseline_full_res_agesexsite, file=paste0(outpath, '/', 'baseline_full_res_agesexsite.txt'), sep = "\t", row.names = FALSE)
write.table(baseline_twins_res_agesexsite, file=paste0(outpath, '/', 'baseline_twins_res_agesexsite.txt'), sep = "\t", row.names = FALSE)
write.table(baseline_full_res_agesexsitepcs, file=paste0(outpath, '/', 'baseline_full_res_agesexsitepcs.txt'), sep = "\t", row.names = FALSE)
write.table(baseline_full_res_agesexsiteeducinc, file=paste0(outpath, '/', 'baseline_full_res_agesexsiteeducinc.txt'), sep = "\t", row.names = FALSE)
write.table(baseline_full_res_agesexsiteeducincpcs, file=paste0(outpath, '/', 'baseline_full_res_agesexsiteeducincpcs.txt'), sep = "\t", row.names = FALSE)
write.table(baseline_twins_res_agesexsiteeducincpcs, file=paste0(outpath, '/', 'baseline_twins_res_agesexsiteeducincpcs.txt'), sep = "\t", row.names = FALSE)