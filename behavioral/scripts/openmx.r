# ------------------------------------------------------------------------------
# Program: mulACEc.R  
#  Author: Hermine Maes
#    Date: 02 25 2016 
#
# Twin Bivariate ACE model to estimate means and (co)variances across multiple groups
# Matrix style model - Raw data - Continuous data
# http://ibg.colorado.edu/cdrom2016/maes/MultivariateAnalysis/
# -------|---------|---------|---------|---------|---------|---------|---------|

rm(list=ls())
# ------------------------------------------------------------------------------
# # Load Libraries & Options
library(OpenMx)
library(reshape)
library(psych)
library(Matrix)
library(stringi)
library('ggplot2')
source("/home/d9smith/projects/random_effects/behavioral/scripts/miFunctions2.R")
source("/home/d9smith/projects/random_effects/behavioral/scripts/twin_functions.R")

mxOption(NULL, "Default optimizer", 'SLSQP') # TODO - what does this do?

# Define the path to tge cmig_tools_utils/r directory
funcpath <- '/home/d9smith/github/cmig_tools_internal/cmig_tools_utils/r'
source(paste0(funcpath, '/', 'loadtxt.R'))

# Load Data
twin_file = "/home/d9smith/projects/random_effects/behavioral/twinfiles/twin_IDs_complete.txt"
pheno_file = "/space/syn50/1/data/ABCD/d9smith/random_effects/behavioral/data/pheno/baseline_twins_res_agesexsite.txt"
# pheno_file = "/space/syn50/1/data/ABCD/d9smith/random_effects/behavioral/data/pheno/twins_res_agesexsiteeducincpcs.txt" 

zyg_file = "/home/d9smith/projects/random_effects/behavioral/twinfiles/ABCD_twins_all.txt"

grm_file = "/home/d9smith/projects/random_effects/behavioral/twinfiles/twins_measured_grm.txt"

twin = loadtxt(twin_file)

grm = read.table(grm_file,header = TRUE,sep=",")

twin <- merge(twin, grm, by=c("IID1","IID2"))
twin_complete = twin[twin$IID1_complete==TRUE & twin$IID2_complete==TRUE,]
twin_complete = twin_complete[,c("IID1","IID2","twin1_genetic_zygosity","measured_grm")]

pheno = loadtxt(pheno_file)
pheno1 = pheno
pheno2 = pheno
names(pheno1)[-(1:2)] = paste0(names(pheno)[-(1:2)],1)
names(pheno2)[-(1:2)] = paste0(names(pheno)[-(1:2)],2)

df = merge(twin_complete, pheno2, by.x=c("IID2"), by.y=c("src_subject_id"))
df = merge(df, pheno1, by.x=c("IID1", "eventname"), by.y=c("src_subject_id", "eventname"))

# Note from DS 2022-08-24:
# The resulting dataframe "df" contains 472 twin pairs, all with complete data for each twin. 
# 271 DZ twin pairs, 201 MZ twin pairs.

# twins = read.table(twin_file, header=T)
# twins_dizyg = read.table(dz_file, sep="\t", header=T)
# twins_monozyg = read.table(mz_file, sep="\t", header=T)
# twin_indiv = unique(c(as.character(twins_dizyg[,'twin1_IID1']), as.character(twins_dizyg[, 'twin2_IID2']),as.character(twins_monozyg[,'twin1_IID1']), as.character(twins_monozyg[, 'twin1_IID2'])))

# This df has 933 twin pairs - DS 2022-08-16
# new df has 676 twin pairs - DS 2022-08-16 (pre-residualized)

# Load nih toolbox tasks
# nihtbx = read.table(nihtbx_file, sep="\t", header=T)
# nihtbx_tsks= c(colnames(nihtbx)[startsWith(colnames(df), 'nihtbx') & endsWith(colnames(df), 'uncorrected')], 'lmt_scr_perc_correct', 'pea_wiscv_trs', 'anthro_height_calc')

# Select relevant individuals
# df = df[twin_indiv,]

# reformat dataframes
# df_untouched = df
# df = df[,]

# Monozygotic coded as 1, Dizygotic coded as 2
df$zyg = NA
df[df$twin1_genetic_zygosity=="Monozygotic",]$zyg = 1
df[df$twin1_genetic_zygosity=="Dizygotic",]$zyg = 2

if (0) {
  # save final list of IDs to pass to FEMA - deprecated as of 2022-08-24
iid1 = df$IID1
iid2 = df$IID2
idlist_for_grm = df[,1:3]
idlist__for_grm_outpath = "/space/syn50/1/data/ABCD/d9smith/random_effects/behavioral/data/pheno/twins_ids_for_grm.txt" 
write.table(idlist_for_grm, file=idlist__for_grm_outpath, sep = "\t", row.names = FALSE)

twins_res_agesex_mx<-rbind(pheno[pheno$src_subject_id %in% iid1,],pheno[pheno$src_subject_id %in% iid2,])
outfile = "/space/syn50/1/data/ABCD/d9smith/random_effects/behavioral/data/pheno/twins_res_agesex_mx.txt" 
write.table(twins_res_agesex_mx, file=outfile, sep = "\t", row.names = FALSE)
}

# Write function for ACE Model
estHerit <- function(data, task, measured_grm=FALSE){
  nt        <- 2                         # number of individuals
  nv        <- 1                         # number of variables
  ntv       <- nv*nt                     # number of total variables
  selVars   <- paste(task,c(rep(1,nv),rep(2,nv)),sep="")

  # Select Data for Analysis
  mzData    <- subset(data, zyg==1, selVars)
  dzData    <- subset(data, zyg==2, selVars)
  nl        <- subset(data, select=c('measured_grm', selVars))

  if (measured_grm == TRUE) {

    # Set Starting Values
    sMu      <- 0                        # start value for means
    sVa      <- .2                       # start value for A
    sVc      <- .2                       # start value for C
    sVe      <- .6                       # start value for E

    # ----------------------------------------------------------------------------------------------------------------------
    # PREPARE MODEL

    # Create Algebra for expected Mean Matrices
    intercept     <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=sMu, labels="interC", name="intercept" )
    expMean       <- mxAlgebra(expression = intercept, name="expMean" )

    # Create Matrices for Variance Components
    covA      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=sVa, label="VA11", name="A" ) 
    covC      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=sVc, label="VC11", name="C" )
    covE      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=sVe, label="VE11", name="E" )

    # Create Algebra for expected Variance/Covariance Matrices in MZ & DZ twins
    covP      <- mxAlgebra( expression= A+C+E, name="V" )
    relA      <- mxMatrix( type="Stand", nrow=nt, ncol=nt, free=FALSE, labels="measured_grm", name="rA" ) 
    relC      <- mxMatrix( type="Unit", nrow=nt, ncol=nt, free=FALSE, name="rC" ) 
    relE      <- mxMatrix( type="Iden", nrow=nt, ncol=nt, free=FALSE, name="rE" ) 

    expCov   <- mxAlgebra( expression= rA%x%A + rC%x%C + rE%x%E, name="expCov" )

    # Create Data Objects for Multiple Groups
    dataTW    <- mxData( observed=nl, type="raw" )

    # Create Expectation Objects for Multiple Groups
    expTW     <- mxExpectationNormal(covariance="expCov", means="intercept", dimnames=selVars )
    funML     <- mxFitFunctionML()

    # Create Model Objects for Multiple Groups
    # defs      <- list(defAge, defSex, relA)
    # pars      <- list( intercept, betaS, betaA, covA, covC, covE, covP, relC, relE)
    defs      <- list(relA)
    pars      <- list(intercept, covA, covC, covE, covP, relC, relE)

    # Create Algebra for Variance Components
    rowVC     <- rep('VC',nv)
    colVC     <- rep(c('A','C','E','SA','SC','SE'),each=nv)
    estVC     <- mxAlgebra( expression=cbind(A,C,E,A/V,C/V,E/V), name="VC", dimnames=list(rowVC,colVC) )

    # Create Confidence Interval Objects
    ciACE     <- mxCI( "VC[1,4:6]" )

    # Build Model with Confidence Intervals
    modelACE   <- mxModel( "oneACEc", defs, pars, expMean, expCov, dataTW, expTW, funML, estVC, ciACE )
} else {
    # Set Starting Values
    svMe      <- rnorm(1)                       # start value for means
    svPa      <- .4                        # start value for path coefficient
    svPe      <- .6                        # start value for path coefficient for e
    lbPa      <- .0001      
    # ACE Model
    # Create Algebra for expected Mean Matrices
    meanG     <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svMe, labels=task, name="meanG" )

    # Create Matrices for Path Coefficients
    pathA     <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=svPa, label="a11", lbound=lbPa, name="a" ) 
    pathC     <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=svPa, label="c11", lbound=lbPa, name="c" )
    pathE     <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=svPe, label="e11", lbound=lbPa, name="e" )

    # Create Algebra for Variance Components
    covA      <- mxAlgebra( expression=a %*% t(a), name="A" )
    covC      <- mxAlgebra( expression=c %*% t(c), name="C" ) 
    covE      <- mxAlgebra( expression=e %*% t(e), name="E" )

    # Create Algebra for expected Variance/Covariance Matrices in MZ & DZ twins
    covP      <- mxAlgebra( expression= A+C+E, name="V" )
    covMZ     <- mxAlgebra( expression= A+C, name="cMZ" )
    covDZ     <- mxAlgebra( expression= 0.5%x%A+ C, name="cDZ" )
    expCovMZ  <- mxAlgebra( expression= rbind( cbind(V, cMZ), cbind(t(cMZ), V)), name="expCovMZ" )
    expCovDZ  <- mxAlgebra( expression= rbind( cbind(V, cDZ), cbind(t(cDZ), V)), name="expCovDZ" )

    # Create Data Objects for Multiple Groups
    dataMZ    <- mxData( observed=mzData, type="raw" )
    dataDZ    <- mxData( observed=dzData, type="raw" )

    # Create Expectation Objects for Multiple Groups
    expMZ     <- mxExpectationNormal( covariance="expCovMZ", means="meanG", dimnames=selVars )
    expDZ     <- mxExpectationNormal( covariance="expCovDZ", means="meanG", dimnames=selVars )
    funML     <- mxFitFunctionML()

    # Create Model Objects for Multiple Groups
    pars      <- list(meanG, pathA, pathC, pathE, covA, covC, covE, covP)
    modelMZ   <- mxModel( name="MZ", pars, covMZ, expCovMZ, dataMZ, expMZ, funML )
    modelDZ   <- mxModel( name="DZ", pars, covDZ, expCovDZ, dataDZ, expDZ, funML )
    multi     <- mxFitFunctionMultigroup( c("MZ","DZ") )

    # Create Algebra for Variance Components
    rowVC     <- rep('VC',nv)
    colVC     <- rep(c('A','C','E','SA','SC','SE', 'SV'),each=nv)
    estVC     <- mxAlgebra( expression=cbind(A,C,E,A/V,C/V,E/V,V), name="VC", dimnames=list(rowVC,colVC))

    # Create Confidence Interval Objects
    ciACE     <- mxCI( "VC[1,1:7]" )

    # Build Model with Confidence Intervals
    modelACE  <- mxModel( "oneACEc", pars, modelMZ, modelDZ, multi, estVC, ciACE )  
  }
  # ------------------------------------------------------------------------------
  # RUN MODEL

  # Run ACE Model
  fitACE    <- mxRun( modelACE, intervals=T )
  sumACE    <- summary( fitACE , verbose=T)

  # Falkoner 
  r_mz = cor(mzData, use='complete.obs')[paste0(task, '1'), paste0(task, '2')]
  r_dz = cor(dzData, use='complete.obs')[paste0(task, '1'), paste0(task, '2')]
#     print(paste0('Exit code status ', fitACE$output$status, ' status status ', fitACE$output$status))
  A <- mxEval(A, fitACE)
  C <- mxEval(C, fitACE)
  E <- mxEval(E, fitACE)

  V <- (A+C+E)   # total variance
  a2 <- A/V      # genetic term as proportion of total variance, i.e. standardized
  c2 <- C/V      # shared environment term as proportion of total variance
  e2 <- E/V      # nonshared environment term as proportion of total variance
  
  # Extract confidence intervals
  CI <- data.frame(
    A=sumACE$CIdetail[sumACE$CIdetail$parameter=='oneACEc.VC[1,4]', 'value'], 
    C=sumACE$CIdetail[sumACE$CIdetail$parameter=='oneACEc.VC[1,5]', 'value'], 
    E=sumACE$CIdetail[sumACE$CIdetail$parameter=='oneACEc.VC[1,6]', 'value'], 
    row.names=c('left', 'right')
  )
  
  c(a2=as.numeric(a2), 
    c2=as.numeric(c2),
    e2=as.numeric(e2),
    CI=CI,
    falkoner=2*(r_mz-r_dz), summary=sumACE)
}

estHerit(df, 'nihtbx_fluidcomp_uncorrected', measured_grm=FALSE)
estHerit(df, 'nihtbx_fluidcomp_uncorrected', measured_grm=TRUE)

estHerit(df, 'nihtbx_cryst_uncorrected', measured_grm=FALSE)

estHerit(df, 'nihtbx_cryst_uncorrected', measured_grm=TRUE)
# Look at warning messages from OpenMx
tasks = names(pheno)[-(1:2)]  

A <- data.frame(
  task=tasks
)
C <- data.frame(
  task=tasks
)
E <- data.frame(
  task=tasks
)
for (t in 1:length(tasks)){
    result = estHerit(df, tasks[t], measured_grm=TRUE)
    A[tasks==tasks[t], 'openmx'] = as.numeric(result$a2)
    A[tasks==tasks[t], 'openmx_ci_lower'] = result$CI.A[1]
    A[tasks==tasks[t], 'openmx_ci_upper'] = result$CI.A[2]
    # C compoenent
    C[tasks==tasks[t], 'openmx'] = as.numeric(result$c2)
    C[tasks==tasks[t], 'openmx_ci_lower'] = result$CI.C[1]
    C[tasks==tasks[t], 'openmx_ci_upper'] = result$CI.C[2]
    # E compoenent
    E[tasks==tasks[t], 'openmx'] = as.numeric(result$e2)
    E[tasks==tasks[t], 'openmx_ci_lower'] = result$CI.E[1]
    E[tasks==tasks[t], 'openmx_ci_upper'] = result$CI.E[2]
}

A$task <- factor(A$task, levels = A$task)
ggplot(A[1:11,]) +
    geom_bar( aes(x=task, y=openmx), stat="identity", fill="skyblue", alpha=0.7)+
    geom_errorbar( aes(x=task, ymin=openmx_ci_lower, ymax=openmx_ci_upper), width=0.4, colour="orange", alpha=0.9, size=1.3) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1), 
         text = element_text(size=20)) + 
    ylab('Heritability')

write.csv(A, '/home/rloughna/data/ABCD/heritability/twins/open_mx/herit_tlbx_residulised_with_CIs.csv',
    row.names=F)

# Uncorrected
ACE = data.frame(task=tasks,E=E$openmx, C=C$openmx, A=A$openmx)
ACE$task <- factor(ACE$task, levels = ACE$task)
ACE = melt(ACE)
ggplot(ACE, aes(x=task, y=value, fill=variable)) +
    geom_bar(stat="identity") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
#     geom_errorbar( aes(x=task, ymin=openmx_ci_lower, ymax=openmx_ci_upper), width=0.4, colour="orange", alpha=0.9, size=1.3) + 

# Residualised for covariates
ACE = data.frame(task=tasks,E=E$openmx, C=C$openmx, A=A$openmx)
ACE$task <- factor(ACE$task, levels = ACE$task)
ACE = melt(ACE)
ggplot(ACE, aes(x=task, y=value, fill=variable)) +
    geom_bar(stat="identity") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
#


