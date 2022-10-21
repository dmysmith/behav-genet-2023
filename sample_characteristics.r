################################
# Create table of demographic information
# Diana Smith
# October 2022

rm(list=ls())
################################
# Create table of sample information
setwd('/home/d9smith/projects/behav_genet_2022')
source('create_phenofiles.R')
library(tableone)

catVars <- c("high.educ", "household.income")

### baseline full sample
t1_full <- CreateTableOne(vars = c("interview_age", catVars), data = baseline_full, factorVars = catVars)

### baseline twin sample
t1_twins <- CreateTableOne(vars = c("interview_age", catVars), data = twinmat, factorVars = catVars)

print(t1_full, quote = TRUE, noSpaces = TRUE)
print(t1_twins, quote = TRUE, noSpaces = TRUE)

# Investigate number of twins in general sample (excluding twin sub-sample)
# Note that create_twin_grm.m contains the code to measure the number of pairs with genetic relatedness > 0.9.
pregnancy_id_file <- '/home/sabad/requests/pregnancy_ID_07182022.csv'
pregnancy_id <- read.table(pregnancy_id_file, sep = ',', header=TRUE)

longitudinal_notwins_orig = longitudinal_notwins
longitudinal_notwins = merge(longitudinal_notwins_orig, pregnancy_id, by.x = 'src_subject_id', by.y = 'pguid')

pregnancy_ids <- unique(longitudinal_notwins$pregnancyID)
counting_twins <- data.frame(pregnancy_ids)
# head(counting_twins)
counting_twins$num_twins = length(unique(pregnancy_id[pregnancy_id$pregnancyID==counting_twins$pregnancy_ids,]$pguid))   

for (row in 1:nrow(counting_twins)) {
    counting_twins[row, 'num_twins'] <- length(unique(longitudinal_notwins[longitudinal_notwins$pregnancyID==counting_twins[row, 'pregnancy_ids'],]$src_subject_id))
}
table(counting_twins$num_twins) 
# counting_twins[counting_twins$num_twins==2,]$pregnancy_ids
# unique(longitudinal_notwins[longitudinal_notwins$pregnancyID==65202,]$src_subject_id)   
