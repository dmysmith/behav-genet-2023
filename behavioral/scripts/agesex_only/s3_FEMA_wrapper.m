%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run FEMA modeling F, A, S, and E on cognition
%% Diana Smith
%% March 2022

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADD ALL ABCD CMIG tools directories to MATLAB path:
addpath(genpath('/home/d9smith/github/cmig_utils'))
addpath(genpath('/home/d9smith/github/FEMA'))
addpath(genpath('/home/d9smith/github/showVol'))
addpath(genpath('/home/d9smith/github/showSurf'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Specify data release

dataRelease = '4.0'; %'3.0' or '4.0'

% run abcdConfig

cfg = abcdConfig('FEMA');
abcd_sync_path=cfg.data.abcd_sync;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Specify where to store results

outDir = '/space/syn50/1/data/ABCD/d9smith/reading/agesex_only'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INPUTS TO FEMA_wrapper.m

dirname_tabulated = fullfile(abcd_sync_path,dataRelease,'tabulated/released'); % directory to tabulated imaging data on abcd-sync 

atlasVersion = 'ABCD2_cor10';
dirname_tabulated = fullfile(abcd_sync_path,'4.0','tabulated/released'); %KNOWN ISSUE: breaks when using txt files following NDA release --> must use pre-release csv files
fname_design = '/home/d9smith/projects/reading/agesex_only/designMat_agesex_only.txt';
fname_pihat = fullfile(abcd_sync_path, '3.0','genomics','ABCD_rel3.0_pihat.mat'); %FIXME: replace with 4.0 when available

% NOTE: `fname_design` can also be a cell array of multiple design matrices
% (txt file paths) that will be looped over in `FEMA_wrapper.m`
% Use makeDesign.R from ~/github/cmig_utils/r to make design matrix compatible with `FEMA_wrapper.m`
% See makeDesign_demo.R for tutorial on how to make design matrices
% compatible with `FEMA_wrapper.m`

% Optional inputs for `FEMA_wrapper.m` depending on analysis
contrasts=[]; % Contrasts relate to columns in design matrix e.g. [1 -1] will take the difference between cols 1 and 2 in your design matrix (X)
ranknorm = 1; % Rank normalizes dependent variables (Y) (default = 0)
nperms = 0; % Number of permutations - if wanting to use resampling methods nperms>0
mediation = 0; % If wanting to use outputs for a mediation analysis set mediation=1 - ensures same resampling scheme used for each model in fname_design
PermType = 'wildbootstrap'; %Default resampling method is null wild-bootstrap - to run mediation analysis need to use non-null wild-bootstrap ('wildboostrap-nn')
tfce_cols = []; % Columns in design matrix to loop over to calculate TFCE - selecting columns of interest improves efficiency - is `isempty(tfce_cols)` FEMA_wrapper will NOT run TFCE

datatype='external'; % can use txt with columns of ANY data type (e.g. ROIs, behavior) - runs mass univaraite LME across every column

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Baseline with random effects F, E

fstem_imaging='baseline_FE';
dirname_imaging = '/space/syn50/1/data/ABCD/d9smith/reading/nih_tbx_baseline.txt';
dirname_out = outDir; %filepath to save FEMA output

RandomEffects = {'F','E'}; % Random effects to include: family, subject, error

% RUN FEMA
[fnames_out beta_hat beta_se zmat logpmat sig2tvec sig2mat beta_hat_perm beta_se_perm zmat_perm sig2tvec_perm sig2mat_perm logLikvec_perm inputs] = FEMA_wrapper(fstem_imaging, fname_design, dirname_out, dirname_tabulated, dirname_imaging, datatype,...
'ranknorm', ranknorm, 'contrasts', contrasts, 'RandomEffects', RandomEffects, 'pihat_file', fname_pihat, 'nperms', nperms, 'mediation',mediation,'PermType',PermType,'tfce_cols',tfce_cols);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Baseline with random effects F, A, E

fstem_imaging='baseline_FAE';
dirname_imaging = '/space/syn50/1/data/ABCD/d9smith/reading/nih_tbx_baseline.txt';
dirname_out = outDir; %filepath to save FEMA output

RandomEffects = {'F','A','E'}; % Random effects to include: family, subject, error

% RUN FEMA
[fnames_out beta_hat beta_se zmat logpmat sig2tvec sig2mat beta_hat_perm beta_se_perm zmat_perm sig2tvec_perm sig2mat_perm logLikvec_perm inputs] = FEMA_wrapper(fstem_imaging, fname_design, dirname_out, dirname_tabulated, dirname_imaging, datatype,...
'ranknorm', ranknorm, 'contrasts', contrasts, 'RandomEffects', RandomEffects, 'pihat_file', fname_pihat, 'nperms', nperms, 'mediation',mediation,'PermType',PermType,'tfce_cols',tfce_cols);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Longitudinal with random effects F, S, E

fstem_imaging='longitudinal_FSE';
dirname_imaging = '/space/syn50/1/data/ABCD/d9smith/reading/nih_tbx_longitudinal.txt';
dirname_out = outDir; %filepath to save FEMA output

RandomEffects = {'F','S','E'}; % Random effects to include: family, subject, error

% RUN FEMA
[fnames_out beta_hat beta_se zmat logpmat sig2tvec sig2mat beta_hat_perm beta_se_perm zmat_perm sig2tvec_perm sig2mat_perm logLikvec_perm inputs] = FEMA_wrapper(fstem_imaging, fname_design, dirname_out, dirname_tabulated, dirname_imaging, datatype,...
'ranknorm', ranknorm, 'contrasts', contrasts, 'RandomEffects', RandomEffects, 'pihat_file', fname_pihat, 'nperms', nperms, 'mediation',mediation,'PermType',PermType,'tfce_cols',tfce_cols);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. Longitudinal with random effects F, A, S, E

fstem_imaging='longitudinal_FASE';
dirname_imaging = '/space/syn50/1/data/ABCD/d9smith/reading/nih_tbx_longitudinal.txt';
dirname_out = outDir; %filepath to save FEMA output

RandomEffects = {'F','A','S','E'}; % Random effects to include: family, subject, error

% RUN FEMA
[fnames_out beta_hat beta_se zmat logpmat sig2tvec sig2mat beta_hat_perm beta_se_perm zmat_perm sig2tvec_perm sig2mat_perm logLikvec_perm inputs] = FEMA_wrapper(fstem_imaging, fname_design, dirname_out, dirname_tabulated, dirname_imaging, datatype,...
'ranknorm', ranknorm, 'contrasts', contrasts, 'RandomEffects', RandomEffects, 'pihat_file', fname_pihat, 'nperms', nperms, 'mediation',mediation,'PermType',PermType,'tfce_cols',tfce_cols);

