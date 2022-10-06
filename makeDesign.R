###############################################################

# Create design matrix to use for FEMA
# Diana Smith
# April 2022

# Note: Because this project involves a pre-residualization step,
# all models will use an "empty" design matrix.
# This means that FEMA will be fitting only the random effects.
###############################################################

# Define path to the makeDesign function (available from:
# https://github.com/cmig-research-group/cmig_tools )
source('/home/d9smith/github/cmig_tools/cmig_tools_utils/r/makeDesign.R')

###############################################################
# The function makeDesign reads variables from an R data frame. This 
# data frame can be loaded as the official ABCD RDS file or any other 
# data frame, for example if you have saved individual instruments as a 
# mini RDS as in the example script createMiniRDS.R

# Load the ABCD data 
ndafile <- '/space/syn50/1/data/ABCD/d9smith/random_effects/behavioral/data/nda4.0_offrel_behavioral.RDS'
nda <- readRDS(ndafile)

# Define the path to the directory where you would like to save out your design matrix 
outpath <- '/space/syn50/1/data/ABCD/d9smith/random_effects/behavioral/designMat'

###############################################################
# Create "empty" design matrix
# load age and sex only file
infile = "/space/syn50/1/data/ABCD/d9smith/random_effects/behavioral/designMat/designMat2_agesex.txt"
agesex = read.table(infile,sep="\t",header=T)

# remove unnecessary rows
agesex = agesex[,c("src_subject_id","eventname","rel_family_id","age","intercept")]

# Define the name of your design matrix file 
fname <- 'designMat0_empty.txt'
# path to output directory
outfile <- paste0(outpath, '/', fname)

# save
write.table(agesex, file=outfile, sep = "\t", row.names = FALSE)