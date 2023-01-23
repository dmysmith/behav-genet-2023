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

# Write a function that takes a dataframe and returns the number of twin/triplet pairs.
count_twins <- function(df, preg_file) {
    preg = read.table(preg_file, sep = ',', header=TRUE)
    df = merge(df, preg, by.x = 'src_subject_id', by.y = 'pguid')
    pregnancy_ids <- unique(df$pregnancyID)
    counting_twins <- data.frame(pregnancy_ids)
    counting_twins$num_twins = length(unique(preg[preg$pregnancyID==counting_twins$pregnancy_ids,]$pguid))   

    for (row in 1:nrow(counting_twins)) {
        counting_twins[row, 'num_twins'] <- length(unique(df[df$pregnancyID==counting_twins[row, 'pregnancy_ids'],]$src_subject_id))
    }

    return(table(counting_twins$num_twins))

}

count_twins(longitudinal_notwins_orig, pregnancy_id_file)

count_twins(baseline_full_res_agesex, pregnancy_id_file)
count_twins(baseline_twins_res_agesex, pregnancy_id_file)
count_twins(longitudinal_full, pregnancy_id_file)


# how many MZ / DZ twins are there in final twin sample?
zygosity_file <- '/home/d9smith/projects/random_effects/behavioral/twinfiles/ABCD_twins_all.txt'
zygosity <- read.table(zygosity_file, header = T, sep = "\t")
twin_unique # df containing the final set of twins
tmp <- merge(twin_unique, zygosity)
table(tmp$twin1_genetic_zygosity)