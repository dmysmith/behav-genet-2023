#calculating heritability of phenotypes based on monozygotic and dizygotic twins at baseline

# libraries
rm(list=ls())
library(data.table)
library(psych)
library(tableone)
library(dplyr)

# define paths
setwd("/home/d9smith/projects/random_effects/")
twinfile = "/home/d9smith/projects/random_effects/behavioral/twinfiles/ABCD_zygozygosity.2020_04_19.check_filt.csv"
phenofile = '/space/syn50/1/data/ABCD/d9smith/random_effects/behavioral/nih_tbx_baseline.txt'
agesexfile = '/home/d9smith/projects/random_effects/behavioral/designMat/designMat2_agesex.txt' 

outpath_mz = "/home/d9smith/projects/random_effects/behavioral/twinfiles/ABCD_mz_IIDs.txt" 
outpath_dz = "/home/d9smith/projects/random_effects/behavioral/twinfiles/ABCD_dz_IIDs.txt"
outpath_all = "/home/d9smith/projects/random_effects/behavioral/twinfiles/ABCD_twins_all.txt" 
outpath_h2 = "/home/d9smith/projects/random_effects/behavioral/twinfiles/twin_h2.txt" 
outpath_twinIDs = "/home/d9smith/projects/random_effects/behavioral/twinfiles/twin_IDs.txt"

pheno<-fread(file=phenofile,stringsAsFactors=FALSE,data.table=FALSE)
zyg<-fread(file=twinfile,stringsAsFactors=FALSE, data.table=FALSE)
phenonames <- names(pheno)

ABCD_twins1<-merge(zyg,pheno,by.x="IID1",by.y="src_subject_id")
colnames(ABCD_twins1) <- paste0('twin1_', colnames(ABCD_twins1))

pheno2 = pheno
colnames(pheno2) <- paste0('twin2_', colnames(pheno))

ABCD_twins<-merge(ABCD_twins1,pheno2,by.x="twin1_IID2",by.y="twin2_src_subject_id")

ABCD_twins <- ABCD_twins %>% rename("IID2" = "twin1_IID2", "IID1" = "twin1_IID1")

twin_ids<-ABCD_twins[,c(1:2)]
ABCD_age<-fread(file=agesexfile,stringsAsFactors=FALSE, data.table=FALSE)
ABCD_age<-ABCD_age[ABCD_age$eventname == "baseline_year_1_arm_1",c("src_subject_id", "age")]
twin_age1<-merge(twin_ids, ABCD_age, by.x=c("IID1"), by.y=c("src_subject_id"))
colnames(twin_age1)[3]<-c("twin1_age")
twin_age2<-merge(twin_ids, ABCD_age, by.x=c("IID2"), by.y=c("src_subject_id"))
colnames(twin_age2)[3]<-c("twin2_age")

twin_age<-merge(twin_age1,twin_age2,by=c("IID1","IID2"))
twin_age$twin<-ifelse(twin_age$twin1_age==twin_age$twin2_age, 1,0)
twin_true<-unique(twin_age[twin_age$twin==1,])

#dz_v2<-dz_true %>% inner_join(dz,by=c("IID1","IID2"))

ABCD_twins<-twin_true %>% inner_join(ABCD_twins,by=c("IID1","IID2"))

#separate by Mono vs dizyg
mz<-ABCD_twins[ABCD_twins$twin1_genetic_zygosity == "Monozygotic",]


#need to remove the dz subjects here that are siblings and not just twins
dz<-ABCD_twins[ABCD_twins$twin1_genetic_zygosity == "Dizygotic",]

write.table(mz, outpath_mz, sep = "\t", quote=FALSE, row.names=FALSE)
write.table(dz, outpath_dz, sep = "\t", quote=FALSE, row.names=FALSE)
write.table(ABCD_twins, outpath_all, sep = "\t", quote=FALSE, row.names=FALSE)

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
  region1<-dz[,c(paste0("twin1_",f))]
  region2<-dz[,c(paste0("twin2_",f))]
  
  corcoef<-cor.test(region1,region2,method="pearson")
  r<-round(corcoef$estimate,3)
  cor[f,2]<-r
  
  #calculate heritability
  h<-2*(cor[f,1]-cor[f,2])
  cor[f,3]<-h
  
}

write.table(cor,outpath_h2, sep = "\t", quote=FALSE, row.names=TRUE) 
#final sample based on 349 monozygotic twins and 491 dizygotic twin pairs

#get final twin sample
twin_IDs<-ABCD_twins[,c(1:2)]

print(paste0("final sample: ", dim(mz)[1], " MZ pairs, ", dim(dz)[1], " DZ pairs"))

# final sample: 349 MZ pairs, 491 DZ pairs

write.table(twin_IDs, outpath_twinIDs, sep = "\t", quote=FALSE, row.names=FALSE) 
