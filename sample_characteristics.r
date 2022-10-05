################################

# Create csv of cognition variables to use for FEMA
# Diana Smith
# April 2022

rm(list=ls())
################################

source('/home/d9smith/projects/random_effects/behavioral/scripts/create_phenofiles.R')
library(tableone)

catVars <- c("high.educ", "household.income")

### baseline full sample
t1_full <- CreateTableOne(vars = c("interview_age", catVars), data = baseline_full, factorVars = catVars)

### baseline twin sample
t1_twins <- CreateTableOne(vars = c("interview_age", catVars), data = twinmat, factorVars = catVars)

print(t1_full, quote = TRUE, noSpaces = TRUE)
print(t1_twins, quote = TRUE, noSpaces = TRUE)
