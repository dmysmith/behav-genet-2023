###############################################################

# Create design matrix of typical covariates to use for FEMA
# Diana Smith
# March 2022

###############################################################
# Define path to the makeDesign function 
source('/home/d9smith/github/cmig_utils/r/makeDesign.R')

###############################################################
# The function makeDesign reads variables from an R data frame. This 
# data frame can be loaded as the official ABCD RDS file or any other 
# data frame, for example if you have saved individual instruments as a 
# mini RDS as in the example script createMiniRDS.R

# Load the data 
ndafile <- '/home/d9smith/FEMA_demo/nda4.0_offrel.RDS'
nda <- readRDS(ndafile)

# Define the path to the directory where you would like to save out your design matrix 
outpath <- '/home/d9smith/projects/reading/agesex_only'

###############################################################
# Design Matrix #1: Cross-sectional

# Define the name of your design matrix file 
fname <- 'designMat_agesex_only.txt'
# path to output directory
outfile <- paste0(outpath, '/', fname)

# makeDesign encodes continuous and categorical variables differently, therefore they are 
# specified using different flags. Continuous variables are added using the "contvar" flag. 
# Here we include age and the genetic PCs as continuous variables.
contvar <- c('interview_age')

# Categorical variables are dummy coded and added using the "catvar" flag. makeDesign automatically
# includes an intercept. For each categorical variable one category is defined as the reference category 
# and that column is dropped. makeDesign also checks whether your matrix is rank deficient and if so
# will automatically drop further categories to avoid this. Here we include sex, scanner info and 
# SES demographics as categorical.
catvar <- c('sex')

# The time points for which we wish to extract data are specified. You do not have to specify multiple 
# time points, however if you do, specify them in chronoligical order ( this will be important for 
# longitudinal modelling)
time <- c('baseline_year_1_arm_1', '2_year_follow_up_y_arm_1') 

# You can specify whether you wish to demean your continuous variables or not. the default is set to 
# demean=TRUE. Other flags include "delta" for longitudinal modelling, "interact" to include interactions, 
# "subjs" if you wish to filter by subject and "quadratic" if you wish to include a quadratic term 
# (more details on these in the following sections). The defaults to these are set to null. 

# Now run makeDesign! 
makeDesign(nda, outfile, time, contvar=contvar, catvar=catvar, delta=NULL, interact=NULL, subjs=NULL, demean=TRUE, quadratic=NULL)
