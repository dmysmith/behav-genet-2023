################################

# Create RDS file for random effects behavioral analysis 
# Diana Smith
# April 2022

################################
# load packages

library(tidyverse)
library(psych)
library(plyr)
library(dplyr)
library(PerformanceAnalytics)
library(pracma)

################################
# Define paths

# tabulated ABCD data 
inpath <- '/space/abcd-sync/4.0/tabulated/released'

# genetic PCs subject data
pcfile <- '/space/gwas-syn2/1/data/GWAS/ABCD/genotype/plink2.eigenvec'

# path to the output RDS file 
outpath <- '/space/syn50/1/data/ABCD/d9smith/random_effects/behavioral'
fname <- 'nda4.0_offrel_behavioral.RDS'
outmatfile <- paste0(outpath, '/', fname)

# cmig_utils/r directory
funcpath <- '/home/d9smith/github/cmig_tools/cmig_tools_utils/r'

# The functionmakeDEAPdemos.R requires the path to the directory which 
# contains the tabulated ABCD data defined explicitly here
deappath <- inpath

# file names for instruments of interest
physfile <- 'abcd_ant01.txt'

# full paths to these files 
physfile <- paste0(inpath, '/', physfile)

################################

# source load.txt and makeDEAPdemos.R 
source(paste0(funcpath, '/', 'loadtxt.R'))
source(paste0(funcpath, '/', 'makeDEAPdemos.R'))

################################
# Load the physical health instrument file  
phys <- loadtxt(physfile)
# Extract the variables of interest
physvar <- c('src_subject_id', 'eventname', 'interview_age', 'sex')
phys <- phys[,physvar]
# Write to a dataframe 
outmat <- phys

################################
# Load the genetic PCs file (this does not require loadtxt.R!) 
pc_mat <- read.delim(pcfile)
# Get just the first 10 PCs and write to a dataframe  
pc <- data.frame(pc_mat[,c('IID','PC1','PC2','PC3','PC4','PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10')])
names(pc) <- c('src_subject_id','PC1','PC2','PC3','PC4','PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10')
# Combine with the physical health variables. 
outmat <- join(outmat, pc, by='src_subject_id', match = "all")

################################
# Create the SES variables as coded by DEAP
deap <- makeDEAPdemos(deappath)
deap <- deap[ , -which(names(deap) %in% c("sex", "interview_date", "interview_age"))]
# Combine with the previously extracted variables
outmat <- join(outmat, deap, by=c('src_subject_id', 'eventname'))

################################
# Save the "outmat" as an RDS 

if ( ! dir.exists(outpath) ) {
        dir.create(outpath, recursive=TRUE)
}


saveRDS(outmat, file=outmatfile)


