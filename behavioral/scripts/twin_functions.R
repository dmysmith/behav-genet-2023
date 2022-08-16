#' Funciton residualise tasks
#'
#' - twins should be a table with columns 'IID1', 'IID2'and 'Zygosity'
#' - df is the abcd nda data release
#' - tasks is a vector with each element being a task of interest (i.e. nihtbx_reading_uncorrected)
residualise_tasks <- function(df, tasks, remove_strs, covariates=c('race.ethnicity',  'sex',  'married',  'interview_age',  'household.income')){
    for (i in 1:length(tasks)){
        # Residualize task
        lmodel <- lm(paste0(tasks[i], ' ~  ', paste(covariates, collapse=' + ')), data=df)
        res = residuals(lmodel)
        # Remove prefix and suffix
        for (remove_str in remove_strs){tasks[i] = gsub(remove_str, "", tasks[i])}
        tasks[i] = paste0(tasks[i], '_residualized')
        df[names(res),tasks[i]] = res
    }
    return (list('df'=df, 'tasks'=tasks))
}



#' Funciton to create nl table, each row is a twin pair with collumns being phenotype1 and phenotype2 e.g. ('reading1' and 'reading2' for each twin)
#'
#' Monozygotic coded as 1, Dizygotic coded as 2
#' - twins should be a table with columns 'IID1', 'IID2'and 'Zygosity'
#' - df is the abcd nda data release
#' - tasks is a vector with each element being a task of interest (i.e. nihtbx_reading_uncorrected)
create_nl <- function(twins, df, tasks){
    df[, tasks] = as.numeric(scale(df[, tasks])) # Z score tasks
    twins_dizyg = twins[twins$Zygosity == 'Dizygotic', ]
    twins_monozyg = twins[twins$Zygosity == 'Monozygotic', ]
    
    nl = data.frame(zyg=c(rep(1, dim(twins_monozyg)[1]), rep(2, dim(twins_dizyg)[1])))
    twin_ind = rbind(twins_monozyg, twins_dizyg)
    for (i in 1:length(tasks)){
      # Add task data to nl - task1 (twin 1) task2 (twin2)
      nl[, paste0(tasks[i], '1')] = df[as.character(twin_ind[,'IID1']), tasks[i]]
      nl[, paste0(tasks[i], '2')] = df[as.character(twin_ind[,'IID2']), tasks[i]]
     # Rescale to have variances around 1
#       variance = var(df[, tasks[i]], use='na.or.complete')
#       nl[, paste0(tasks[i], '1')] = df[as.character(twin_ind[,'IID1']), tasks[i]]/variance^.5
#       nl[, paste0(tasks[i], '2')] = df[as.character(twin_ind[,'IID2']), tasks[i]]/variance^.5
    }
    # Omit NAn Entries
    nl= nl[complete.cases(nl), ]
}




#' Function to create ACE model
#'
#' @param nl dataframe Phenotype matrix for matched twins
#' @param var string The name of the phenotype
create_sat_model <- function(nl, vars){
    # Select Variables for Analysis
    nv        <- 1                         # number of variables
    ntv       <- nv*2                      # number of total variables
    selVars   <- paste(vars,c(rep(1,nv),rep(2,nv)),sep="")

    # Select Data for Analysis
    mzData    <- subset(nl, zyg==1, selVars)
    dzData    <- subset(nl, zyg==2, selVars)

    # Set Starting Values
    svMe      <- 5                        # start value for means
    svVa      <- .8                        # start value for variance
    svVas     <- diag(svVa,ntv,ntv)        # assign start values to diagonal of matrix
    lbVa      <- .0001                     # start value for lower bounds
    lbVas     <- diag(lbVa,ntv,ntv)        # assign lower bounds values to diagonal of matrix

    # ------------------------------------------------------------------------------
    # PREPARE MODEL

    # Saturated Model
    # Create Algebra for expected Mean Matrices
    meanMZ    <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svMe, labels=c("mMZ1","mMZ2"), name="meanMZ" )
    meanDZ    <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svMe, labels=c("mDZ1","mDZ2"), name="meanDZ" )

    # Create Algebra for expected Variance/Covariance Matrices
    covMZ     <- mxMatrix( type="Symm", nrow=ntv, ncol=ntv, free=TRUE, values=svVas, lbound=lbVas, labels=c("vMZ1","cMZ21","vMZ2"), name="covMZ" )
    covDZ     <- mxMatrix( type="Symm", nrow=ntv, ncol=ntv, free=TRUE, values=svVas, lbound=lbVas, labels=c("vDZ1","cDZ21","vDZ2"), name="covDZ" )

    # Create Algebra for Maximum Likelihood Estimates of Twin Correlations
    matI      <- mxMatrix( type="Iden", nrow=ntv, ncol=ntv, name="I" )
    corMZ     <- mxAlgebra( solve(sqrt(I*covMZ)) %&% covMZ, name="corMZ" )
    corDZ     <- mxAlgebra( solve(sqrt(I*covDZ)) %&% covDZ, name="corDZ" )

    # Create Data Objects for Multiple Groups
    dataMZ    <- mxData( observed=mzData, type="raw" )
    dataDZ    <- mxData( observed=dzData, type="raw" )

    # Create Expectation Objects for Multiple Groups
    expMZ     <- mxExpectationNormal( covariance="covMZ", means="meanMZ", dimnames=selVars )
    expDZ     <- mxExpectationNormal( covariance="covDZ", means="meanDZ", dimnames=selVars )
    funML     <- mxFitFunctionML()

    # Create Model Objects for Multiple Groups
    modelMZ   <- mxModel( name="MZ", meanMZ, covMZ, matI, corMZ, dataMZ, expMZ, funML )
    modelDZ   <- mxModel( name="DZ", meanDZ, covDZ, matI, corDZ, dataDZ, expDZ, funML )
    multi     <- mxFitFunctionMultigroup( c("MZ","DZ") )

    # Create Confidence Interval Objects
    ciCov     <- mxCI( c('MZ.covMZ','DZ.covDZ') )
    ciMean    <- mxCI( c('MZ.meanMZ','DZ.meanDZ') )

    # Build Saturated Model with Confidence Intervals
    model     <- mxModel( "oneSATc", modelMZ, modelDZ, multi, ciCov, ciMean )
}




#' Function to create ACE model
#'
#' @param nl dataframe Phenotype matrix for matched twins
#' @param var string The name of the phenotype
create_univACE_model <- function(nl, vars){
    # Select Variables for Analysis
    nv        <- 1                         # number of variables
    ntv       <- nv*2                      # number of total variables
    selVars   <- paste(vars,c(rep(1,nv),rep(2,nv)),sep="")

    # Select Data for Analysis
    mzData    <- subset(nl, zyg==1, selVars)
    dzData    <- subset(nl, zyg==2, selVars)

    # Set Starting Values
    svMe      <- rnorm(1)                       # start value for means
    svPa      <- .4                        # start value for path coefficient
    svPe      <- .6                        # start value for path coefficient for e
    lbPa      <- .0001                     # start value for lower bounds

    # ------------------------------------------------------------------------------
    # PREPARE MODEL

    # ACE Model
    # Create Algebra for expected Mean Matrices
    meanG     <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svMe, labels=vars, name="meanG" )

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
    colVC     <- rep(c('A','C','E','SA','SC','SE','V'),each=nv)
    estVC     <- mxAlgebra( expression=cbind(A,C,E,A/V,C/V,E/V,V), name="VC", dimnames=list(rowVC,colVC))

    # Create Confidence Interval Objects
    ciACE     <- mxCI( "VC[1,1:7]" )

    # Build Model with Confidence Intervals
    modelACE  <- mxModel( "oneACEc", pars, modelMZ, modelDZ, multi, estVC, ciACE )
    return(modelACE)
}


#' Function to create ADE model
#'
#' @param nl dataframe Phenotype matrix for matched twins
#' @param var string The name of the phenotype
create_univADE_model <- function(nl, vars){
    # Select Variables for Analysis
    nv        <- 1                         # number of variables
    ntv       <- nv*2                      # number of total variables
    selVars   <- paste(vars,c(rep(1,nv),rep(2,nv)),sep="")

    # Select Data for Analysis
    mzData    <- subset(nl, zyg==1, selVars)
    dzData    <- subset(nl, zyg==2, selVars)

    # Set Starting Values
    svMe      <- rnorm(1)                       # start value for means
    svPa      <- .4                        # start value for path coefficient
    svPe      <- .6                        # start value for path coefficient for e
    lbPa      <- .0001                     # start value for lower bounds
    # ------------------------------------------------------------------------------
    # PREPARE MODEL

    # ADE Model
    # Create Algebra for expected Mean Matrices
    meanG     <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svMe, labels="xbmi", name="meanG" )

    # Create Matrices for Path Coefficients
    pathA     <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=svPa, label="a11", lbound=lbPa, name="a" ) 
    pathD     <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=svPa, label="d11", lbound=lbPa, name="d" )
    pathE     <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=svPe, label="e11", lbound=lbPa, name="e" )

    # Create Algebra for Variance Components
    covA      <- mxAlgebra( expression=a %*% t(a), name="A" )
    covD      <- mxAlgebra( expression=d %*% t(d), name="D" ) 
    covE      <- mxAlgebra( expression=e %*% t(e), name="E" )

    # Create Algebra for expected Variance/Covariance Matrices in MZ & DZ twins
    covP      <- mxAlgebra( expression= A+D+E, name="V" )
    covMZ     <- mxAlgebra( expression= A+D, name="cMZ" )
    covDZ     <- mxAlgebra( expression= 0.5%x%A+ 0.25%x%D, name="cDZ" )
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
    pars      <- list(meanG, pathA, pathD, pathE, covA, covD, covE, covP)
    modelMZ   <- mxModel( name="MZ", pars, covMZ, expCovMZ, dataMZ, expMZ, funML )
    modelDZ   <- mxModel( name="DZ", pars, covDZ, expCovDZ, dataDZ, expDZ, funML )
    multi     <- mxFitFunctionMultigroup( c("MZ","DZ") )

    # Create Algebra for Variance Components
    rowVC     <- rep('VC',nv)
    colVC     <- rep(c('A','D','E','SA','SD','SE'),each=nv)
    estVC     <- mxAlgebra( expression=cbind(A,D,E,A/V,D/V,E/V), name="VC", dimnames=list(rowVC,colVC))

    # Create Confidence Interval Objects
    ciADE     <- mxCI( "VC[1,1:3]" )

    # Build Model with Confidence Intervals
    modelADE  <- mxModel( "oneADEc", pars, modelMZ, modelDZ, multi, estVC, ciADE )
}

