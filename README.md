# Heritability estimation of cognitive phenotypes in the ABCD Study using mixed models
Diana M. Smith, Rob Loughnan, Naomi Friedman, Wes Thompson, Ole Andreassen, Michael Neale, Terry L. Jernigan, Anders M. Dale.

(Under Review.)

This repository contains the code used to run all analyses in the manuscript. A description of each file is below, listed in the recommended order for running the code.

## create_phenofiles.r
Collects the ABCD data into a set of dataframes to be used as the dependent variables in all models. Several files are saved out (because some models use different samples, adjust for all covariates versus just age and sex, etc).
