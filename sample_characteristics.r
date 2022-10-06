################################
# Create table of demographic information
# Diana Smith
# October 2022

rm(list=ls())
################################

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
