# Heritability estimation of cognitive phenotypes in the ABCD Study using mixed models
Diana M. Smith, Rob Loughnan, Naomi Friedman, Wes Thompson, Ole Andreassen, Michael Neale, Terry L. Jernigan, Anders M. Dale.

(Under Review.)

This repository contains the code used to run all analyses in the manuscript. A description of each file is below, listed in the recommended order for running the code.

## create_phenofiles.r
Collects the ABCD data into a set of dataframes to be used as the dependent variables in all models. Several files are saved out (because some models use different samples, adjust for all covariates versus just age and sex, etc).

## openmx.r
Runs the ACE model in OpenMx and saves the estimates for A, C, and E.

## sample_characteristics.r
Creates table of demographic information.

## makeDesign.R
Creates design matrix to pass to FEMA. For this manuscript the design matrix is "empty," i.e. contains no fixed effects.

## create_twin_grm.m
Creates matrix of "assigned" genetic relatedness values.

## bg_fema_wrapper.m
Script that creates several models to run in FEMA then calls FEMA helper functions (available at https://github.com/cmig-research-group/cmig_tools)

# twin_functions.R and miFunctions2.R
Helper functions for OpenMx analysis.

# visualize_python.ipynb
Code for creating the plots shown in the manuscript.