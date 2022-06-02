################################

# Create csv of cognition variables to use for FEMA
# Diana Smith
# April 2022

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
inpath <- '/space/abcd-sync/4.0/tabulated/released'

# Define the full path to the output RDS file 
outpath <- '/space/syn50/1/data/ABCD/d9smith/random_effects/behavioral/agesex_only'

# Define the path to tge cmig_tools_utils/r directory
funcpath <- '/home/d9smith/github/cmig_tools/cmig_tools_utils/r'
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