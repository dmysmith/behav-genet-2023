# create design matrix for baseline twin analysis

designmat_path = "/home/d9smith/projects/random_effects/behavioral/designMat2_agesex.txt"
outfile = "/home/d9smith/projects/random_effects/behavioral/designMat5_agesextwinsbaseline.txt" 

designmat = fread(file=designmat_path,stringsAsFactors=FALSE,data.table=FALSE)

twin_ids = fread(file="/home/d9smith/projects/random_effects/twinfiles/twin_IDs.txt",stringsAsFactors=FALSE,data.table=FALSE)
twin_id_list = c(twin_ids[,1], twin_ids[,2])

designmat_baseline = designmat[designmat$eventname == "baseline_year_1_arm_1",]
designmat_baseline_twins = designmat_baseline[designmat_baseline$src_subject_id %in% twin_id_list,]

# write file
write.table(designmat_baseline_twins, outfile, col.names=TRUE, row.names=FALSE, sep='\t') 
