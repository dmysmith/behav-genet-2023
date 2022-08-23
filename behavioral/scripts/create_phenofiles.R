################################

# Create csv of cognition variables to use for FEMA
# Diana Smith
# April 2022

rm(list=ls())
################################
# The following R packages need to be loaded

library(tidyverse)
library(psych)
library(plyr)
library(dplyr)
library(PerformanceAnalytics)
library(pracma)

################################
# This section defines input and output paths for files and functions called. 

# Define the path to the directory which contains the tabulated ABCD data 
inpath <- '/space/amdale/1/tmp/ABCD_cache/abcd-sync/4.0/tabulated/released'

# Define the full path to the output RDS file 
outpath <- '/space/syn50/1/data/ABCD/d9smith/random_effects/behavioral/data/pheno'

# Define the path to tge cmig_tools_utils/r directory
funcpath <- '/home/d9smith/github/cmig_tools_internal/cmig_tools_utils/r'
# The functionmakeDEAPdemos.R requires the path to the directory which 
# contains the tabulated ABCD data defined explicitly here
deappath <- inpath

# Define the file names for the instruments from which we need to pull
# variables. 
tbx_file <- 'abcd_tbss01.txt'
lmt_file <- 'lmtp201.txt'
reasoning_file <- 'abcd_ps01.txt'
phys_file <- 'abcd_ant01.txt'
# Define the full paths to these files 
tbx_file <- paste0(inpath, '/', tbx_file)
lmt_file <- paste0(inpath, '/', lmt_file)
reasoning_file <- paste0(inpath, '/', reasoning_file)
phys_file <- paste0(inpath, '/', phys_file)

################################

# R needs to parse two functions from cmig_utils/r. One is load.txt
# to load the .txt files and remove unnecessary 2nd row, the other is 
# makeDEAPdemos.R which collates SES variables into the format used by
# DESP. 
source(paste0(funcpath, '/', 'loadtxt.R'))
source(paste0(funcpath, '/', 'makeDEAPdemos.R'))

################################
# Load the NIH toolbox instrument file  
tbx <- loadtxt(tbx_file)
# Extract the variables of interest
tbxvar <- c('src_subject_id', 'eventname', 'nihtbx_picvocab_uncorrected', 'nihtbx_flanker_uncorrected', 
'nihtbx_list_uncorrected', 'nihtbx_cardsort_uncorrected','nihtbx_pattern_uncorrected','nihtbx_picture_uncorrected', 
'nihtbx_reading_uncorrected', 'nihtbx_fluidcomp_uncorrected','nihtbx_cryst_uncorrected', 'nihtbx_totalcomp_uncorrected')
tbx <- tbx[,tbxvar]
# Write to a dataframe 
outmat <- tbx

################################
# Load the little man task file
lmt <- loadtxt(lmt_file)
# Extract the variables of interest
lmtvar <- c('src_subject_id', 'eventname', 'lmt_scr_perc_correct')
lmt <- lmt[,lmtvar]
# Combine with the physical health variables. 
outmat <- join(outmat, lmt, by=c('src_subject_id', 'eventname'))

################################
# Load the reasoning file
reasoning <- loadtxt(reasoning_file)
# Extract the variables of interest
reasoningvar <- c('src_subject_id', 'eventname', 'pea_wiscv_trs')
reasoning <- reasoning[,reasoningvar]
# Combine with the physical health variables. 
outmat <- join(outmat, reasoning, by=c('src_subject_id', 'eventname'))

################################
# Load the physical variables file  
phys <- loadtxt(phys_file)
# Extract the variables of interest
physvar <- c('src_subject_id', 'eventname', 'anthroheightcalc')
phys <- phys[,physvar]
# Combine with outmat. 
outmat <- join(outmat, phys, by=c('src_subject_id', 'eventname'))

################################
# Create dataframe "baseline" that includes all variables at baseline
baseline <- outmat[outmat$eventname=='baseline_year_1_arm_1',]
baselinevars <- names(outmat)[-(1:2)]
################################
# Create dataframe "longitudinal" that includes baseline and year 2 for all variables with data
y2vars <- c('src_subject_id','eventname', 'nihtbx_picvocab_uncorrected','nihtbx_flanker_uncorrected','nihtbx_pattern_uncorrected',
'nihtbx_picture_uncorrected','nihtbx_reading_uncorrected','nihtbx_cryst_uncorrected','lmt_scr_perc_correct','anthroheightcalc')

longitudinal <- outmat[,y2vars]

################################
# Save both "baseline" and "longitudinal" as RDS files 

if ( ! dir.exists(outpath) ) {
        dir.create(outpath, recursive=TRUE)
}

write.table(baseline, file=paste0(outpath, '/', 'nih_tbx_baseline.txt'), sep = "\t", row.names = FALSE)
write.table(longitudinal, file=paste0(outpath, '/', 'nih_tbx_longitudinal.txt'), sep = "\t", row.names = FALSE)

################################
# Create several residualized .csv files
################################
# Create the SES variables as coded by DEAP
deap <- makeDEAPdemos(deappath)
deap <- deap[ , c("src_subject_id", "eventname", "sex", "interview_age", "high.educ", "household.income")]
# Combine with the previously extracted variables
fullmat <- join(outmat, deap, by=c('src_subject_id', 'eventname'))

# genetic PCs subject data
pcfile <- '/space/gwas-syn2/1/data/GWAS/ABCD/genotype_proc/imputation/pop_struct_smokescreen/ABCD_20220428.updated.nodups.curated_pcair.tsv'
pc_mat <- read.delim(pcfile)
# Get just the first 10 PCs and write to a dataframe  
pc <- data.frame(pc_mat[,c('X','C1','C2','C3','C4','C5', 'C6', 'C7', 'C8', 'C9', 'C10')])
names(pc) <- c('src_subject_id','PC1','PC2','PC3','PC4','PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10')
# Combine with the physical health variables. 
fullmat <- join(fullmat, pc, by='src_subject_id', match = "all")
fullmat_untouched = fullmat
fullmat = na.omit(fullmat)

# 1. all_res_agesex
all_res_agesex = fullmat[,c("src_subject_id", "eventname", baselinevars)]
allModelsList <- lapply(paste(baselinevars, "~ interview_age + sex"), as.formula)
allModelsResults <- lapply(allModelsList, function(x) lm(x, data = fullmat, na.action = na.exclude))
allModelsResiduals <- lapply(allModelsList, function(x) residuals(lm(x, data = fullmat)))  
all_res_agesex[,-(1:2)] = allModelsResiduals

# 2. twins_res_agesex
twin_ids = fread(file="/home/d9smith/projects/random_effects/behavioral/twinfiles/twin_IDs.txt",stringsAsFactors=FALSE,data.table=FALSE)
twin_id_list = c(twin_ids[,1], twin_ids[,2])
twinmat = fullmat[fullmat$src_subject_id %in% twin_id_list,]
twins_res_agesex = twinmat[, c("src_subject_id", "eventname", baselinevars)]

allModelsList <- lapply(paste(baselinevars, "~ interview_age + sex"), as.formula)
allModelsResults <- lapply(allModelsList, function(x) lm(x, data = twinmat, na.action = na.exclude))
allModelsResiduals <- lapply(allModelsList, function(x) residuals(lm(x, data = twinmat)))  
twins_res_agesex[,-(1:2)] = allModelsResiduals

# 3. all_res_agesexpcs
all_res_agesexpcs = fullmat[,c("src_subject_id", "eventname", baselinevars)]
allModelsList <- lapply(paste(baselinevars, "~ interview_age + sex + PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10"), as.formula)
allModelsResults <- lapply(allModelsList, function(x) lm(x, data = fullmat, na.action = na.exclude))
allModelsResiduals <- lapply(allModelsList, function(x) residuals(lm(x, data = fullmat)))  
all_res_agesexpcs[,-(1:2)] = allModelsResiduals

# 4. all_res_agesexeducinc
all_res_agesexeducinc = fullmat[,c("src_subject_id", "eventname", baselinevars)]
allModelsList <- lapply(paste(baselinevars, "~ interview_age + sex + high.educ + household.income"), as.formula)
allModelsResults <- lapply(allModelsList, function(x) lm(x, data = fullmat, na.action = na.exclude))
allModelsResiduals <- lapply(allModelsList, function(x) residuals(lm(x, data = fullmat)))  
all_res_agesexeducinc[,-(1:2)] = allModelsResiduals

# 5. all_res_agesexeducincpcs
all_res_agesexeducincpcs = fullmat[,c("src_subject_id", "eventname", baselinevars)]
allModelsList <- lapply(paste(baselinevars, "~ interview_age + sex + high.educ + household.income + PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10"), as.formula)
allModelsResults <- lapply(allModelsList, function(x) lm(x, data = fullmat, na.action = na.exclude))
allModelsResiduals <- lapply(allModelsList, function(x) residuals(lm(x, data = fullmat)))  
all_res_agesexeducincpcs[,-(1:2)] = allModelsResiduals

# 6. twins_res_agesexeducincpcs
twins_res_agesexeducincpcs = twinmat[,c("src_subject_id", "eventname", baselinevars)]
allModelsList <- lapply(paste(baselinevars, "~ interview_age + sex + high.educ + household.income + PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10"), as.formula)
allModelsResults <- lapply(allModelsList, function(x) lm(x, data = twinmat, na.action = na.exclude))
allModelsResiduals <- lapply(allModelsList, function(x) residuals(lm(x, data = twinmat)))  
twins_res_agesexeducincpcs[,-(1:2)] = allModelsResiduals 

# save phenofiles - all baseline for now - DS 2022-08-16
write.table(all_res_agesex, file=paste0(outpath, '/', 'all_res_agesex.txt'), sep = "\t", row.names = FALSE)
write.table(twins_res_agesex, file=paste0(outpath, '/', 'twins_res_agesex.txt'), sep = "\t", row.names = FALSE)
write.table(all_res_agesexpcs, file=paste0(outpath, '/', 'all_res_agesexpcs.txt'), sep = "\t", row.names = FALSE)
write.table(all_res_agesexeducinc, file=paste0(outpath, '/', 'all_res_agesexeducinc.txt'), sep = "\t", row.names = FALSE)
write.table(all_res_agesexeducincpcs, file=paste0(outpath, '/', 'all_res_agesexeducincpcs.txt'), sep = "\t", row.names = FALSE)
write.table(twins_res_agesexeducincpcs, file=paste0(outpath, '/', 'twins_res_agesexeducincpcs.txt'), sep = "\t", row.names = FALSE)