#calculating heritability of phenotypes based on monozygotic and dizygotic twins at baseline

# libraries
rm(list=ls())
library(data.table)
library(psych)
library(tableone)
library(dplyr)

# define paths
setwd("/home/d9smith/projects/random_effects/")
twinfile = "/home/d9smith/projects/random_effects/twinfiles/from_carolina/ABCD_zygozygosity.2020_04_19.check_filt.csv"
phenofile = '/space/syn50/1/data/ABCD/d9smith/random_effects/behavioral/nih_tbx_baseline.txt'
agesexfile = '/home/d9smith/projects/random_effects/behavioral/designMat2_agesex.txt' 

outpath_mz = "/home/d9smith/projects/random_effects/twinfiles/ABCD_mz_IIDs.txt" 
outpath_dz = "/home/d9smith/projects/random_effects/twinfiles/ABCD_dz_IIDs.txt"
outpath_h2 = "/home/d9smith/projects/random_effects/twinfiles/twin_h2.txt" 
outpath_twinIDs = "/home/d9smith/projects/random_effects/twinfiles/twin_IDs.txt"

pheno<-fread(file=phenofile,stringsAsFactors=FALSE,data.table=FALSE)
zyg<-fread(file=twinfile,stringsAsFactors=FALSE, data.table=FALSE)
phenonames <- names(pheno)

ABCD_twins1<-merge(zyg,pheno,by.x="IID1",by.y="src_subject_id")
colnames(ABCD_twins1) <- paste0('twin1_', colnames(ABCD_twins1))

ABCD_twins2<-merge(zyg,pheno,by.x="IID2",by.y="src_subject_id")
colnames(ABCD_twins2) <- paste0('twin2_', colnames(ABCD_twins2))

#keep only IDs that have phenotype data for both IID1 and IID2 of the same family
ABCD_twins<-merge(ABCD_twins1,ABCD_twins2,by.x="twin1_IID1",by.y="twin2_IID1")

#separate by Mono vs dyzog
mz<-ABCD_twins[ABCD_twins$twin1_genetic_zygosity == "Monozygotic",]


#need to remove the dz subjects here that are siblings and not just twins
dz<-ABCD_twins[ABCD_twins$twin1_genetic_zygosity == "Dizygotic",]

dz_ids<-dz[,c(1:2)]
ABCD_age<-fread(file=agesexfile,stringsAsFactors=FALSE, data.table=FALSE)
ABCD_age<-ABCD_age[ABCD_age$eventname == "baseline_year_1_arm_1",c("src_subject_id", "age")]
dz_age1<-merge(dz_ids, ABCD_age, by.x=c("twin1_IID1"), by.y=c("src_subject_id"))
colnames(dz_age1)[3]<-c("twin1_age")

dz_age2<-merge(dz_ids, ABCD_age, by.x=c("twin1_IID2"), by.y=c("src_subject_id"))
colnames(dz_age2)[3]<-c("twin2_age")

dz_age<-merge(dz_age1,dz_age2,by=c("twin1_IID1","twin1_IID2"))
dz_age$twin<-ifelse(dz_age$twin1_age==dz_age$twin2_age, 1,0)
dz_true<-unique(dz_age[dz_age$twin==1,])

dz_v2<-dz_true %>% inner_join(dz,by=c("twin1_IID1","twin1_IID2"))

write.table(mz, outpath_mz, sep = "\t", quote=FALSE, row.names=FALSE)
write.table(dz_v2, outpath_dz, sep = "\t", quote=FALSE, row.names=FALSE)

#calculate correlation for each Feature in mono vs dizygotic
cor<- data.frame(matrix(ncol = 3, nrow = length(phenonames[-2:0])))
row.names(cor)<-phenonames[-2:0]  
x <- c("r_mz","r_dz","h2")
colnames(cor) <- x

for (f in phenonames[-2:0]){
  
  #monozygotic
  pheno1<-mz[,c(paste0("twin1_",f))]
  pheno2<-mz[,c(paste0("twin2_",f))]
  
  corcoef<-cor.test(pheno1,pheno2,method="pearson")
  r<-round(corcoef$estimate,3)
  cor[f,1]<-r
  
  #dizygotic
  region1<-dz_v2[,c(paste0("twin1_",f))]
  region2<-dz_v2[,c(paste0("twin2_",f))]
  
  corcoef<-cor.test(region1,region2,method="pearson")
  r<-round(corcoef$estimate,3)
  cor[f,2]<-r
  
  #calculate heritability
  h<-2*(cor[f,1]-cor[f,2])
  cor[f,3]<-h
  
}

write.table(cor,outpath_h2, sep = "\t", quote=FALSE, row.names=TRUE) 
#final sample based on 248 monozygotic twins and 374 dizygotic twin pairs

#get final twin sample
twin_IDs<-ABCD_twins[,c(1:2)]
names(twin_IDs)<-c("IID1","IID2")

print(paste0("final sample: ", dim(mz)[1], " MZ pairs, ", dim(dz_v2)[1], " DZ pairs"))

# final sample: 355 MZ pairs, 528 DZ pairs

write.table(twin_IDs, outpath_twinIDs, sep = "\t", quote=FALSE, row.names=FALSE) 
