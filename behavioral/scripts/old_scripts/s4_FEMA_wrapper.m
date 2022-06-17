%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run FEMA modeling F, A, S, and E on cognition
%% Diana Smith
%% March 2022

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADD CMIG tools directory to MATLAB path:
addpath(genpath('/home/d9smith/github/cmig_tools'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Specify data release

dataRelease = '4.0'; %'3.0' or '4.0'

% run abcdConfig

% cfg = abcdConfig('FEMA');
% abcd_sync_path=cfg.data.abcd_sync;
abcd_sync_path='/space/abcd-sync/';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Specify where to store results

outDir = '/space/syn50/1/data/ABCD/d9smith/random_effects/behavioral'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INPUTS TO FEMA_wrapper.m

% dirname_tabulated = fullfile(abcd_sync_path,dataRelease,'tabulated/released'); % directory to tabulated imaging data on abcd-sync 

atlasVersion = 'ABCD2_cor10';
dirname_tabulated = fullfile('/space/amdale/1/tmp/ABCD_cache/abcd-sync/4.0/tabulated/released/'); %KNOWN ISSUE: breaks when using txt files following NDA release --> must use pre-release csv files
%To run multiple deisgn matrices with same imaging data populate each row with path to each design matrix
designmat_dir = '/home/d9smith/projects/random_effects/behavioral';
designmat_file = dir(sprintf('%s/designMat*.txt', designmat_dir));
designmat_file = {designmat_file.name}';
fname_design = strcat(designmat_dir, '/', designmat_file);
% fname_pihat = fullfile(abcd_sync_path, '3.0','genomics','ABCD_rel3.0_pihat.mat'); %FIXME: replace with 4.0 when available
fname_pihat = fullfile('/space/amdale/1/tmp/ABCD_cache/abcd-sync/3.0/genomics/ABCD_rel3.0_pihat.mat');

% NOTE: `fname_design` can also be a cell array of multiple design matrices
% (txt file paths) that will be looped over in `FEMA_wrapper.m`
% Use makeDesign.R from ~/github/cmig_utils/r to make design matrix compatible with `FEMA_wrapper.m`
% See makeDesign_demo.R for tutorial on how to make design matrices
% compatible with `FEMA_wrapper.m`

% Optional inputs for `FEMA_wrapper.m` depending on analysis
contrasts=[]; % Contrasts relate to columns in design matrix e.g. [1 -1] will take the difference between cols 1 and 2 in your design matrix (X)
ranknorm = 0; % Rank normalizes dependent variables (Y) (default = 0)
nperms = 1000; % Number of permutations - if wanting to use resampling methods nperms>0 and load Anders' script for FEMA_fit
mediation = 0; % If wanting to use outputs for a mediation analysis set mediation=1 - ensures same resampling scheme used for each model in fname_design
PermType = 'wildbootstrap'; %Default resampling method is null wild-bootstrap - to run mediation analysis need to use non-null wild-bootstrap ('wildboostrap-nn')
tfce = 0; % Columns in design matrix to loop over to calculate TFCE - selecting columns of interest improves efficiency - is `isempty(tfce_cols)` FEMA_wrapper will NOT run TFCE

datatype='external'; % can use txt with columns of ANY data type (e.g. ROIs, behavior) - runs mass univaraite LME across every column

% if nperms > 0 will need to add Anders' script
if nperms > 0
    addpath(genpath('~dale/matlab/SSE_local'));
    outDir = sprintf('%s/%.0fperms', outDir, nperms);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Baseline with random effects F, E -- not running due to error with <3 random effects
if 0
    fstem_imaging='baseline_FE';
    dirname_imaging = '/space/syn50/1/data/ABCD/d9smith/random_effects/behavioral/nih_tbx_baseline.txt';
    % dirname_out = outDir; %filepath to save FEMA output
    out_path = strcat(outDir, '/', designmat_file);
    dirname_out = out_path;
    
    RandomEffects = {'F','E'}; % Random effects to include: family, subject, error
    
    % RUN FEMA
    [fnames_out beta_hat beta_se zmat logpmat sig2tvec sig2mat beta_hat_perm beta_se_perm zmat_perm sig2tvec_perm sig2mat_perm logLikvec_perm inputs] = FEMA_wrapper(fstem_imaging, fname_design, dirname_out, dirname_tabulated, dirname_imaging, datatype,...
    'ranknorm', ranknorm, 'contrasts', contrasts, 'RandomEffects', RandomEffects, 'pihat_file', fname_pihat, 'nperms', nperms, 'mediation',mediation,'PermType',PermType,'tfce',tfce);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Baseline with random effects F, A, E

fstem_imaging='baseline_FAE';
dirname_imaging = '/space/syn50/1/data/ABCD/d9smith/random_effects/behavioral/nih_tbx_baseline.txt';
% dirname_out = outDir; %filepath to save FEMA output

out_path = strcat(outDir, '/', designmat_file);
dirname_out = out_path;

RandomEffects = {'F','A','E'}; % Random effects to include: family, subject, error

% RUN FEMA
[fnames_out beta_hat beta_se zmat logpmat sig2tvec sig2mat beta_hat_perm beta_se_perm zmat_perm sig2tvec_perm sig2mat_perm logLikvec_perm inputs] = FEMA_wrapper(fstem_imaging, fname_design, dirname_out, dirname_tabulated, dirname_imaging, datatype,...
'ranknorm', ranknorm, 'contrasts', contrasts, 'RandomEffects', RandomEffects, 'pihat_file', fname_pihat, 'nperms', nperms, 'mediation',mediation,'PermType',PermType,'tfce',tfce);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Longitudinal with random effects F, S, E

fstem_imaging='longitudinal_FSE';
dirname_imaging = '/space/syn50/1/data/ABCD/d9smith/random_effects/behavioral/nih_tbx_longitudinal.txt';
% dirname_out = outDir; %filepath to save FEMA output

out_path = strcat(outDir, '/', designmat_file);
dirname_out = out_path;

RandomEffects = {'F','S','E'}; % Random effects to include: family, subject, error

% RUN FEMA
[fnames_out beta_hat beta_se zmat logpmat sig2tvec sig2mat beta_hat_perm beta_se_perm zmat_perm sig2tvec_perm sig2mat_perm logLikvec_perm inputs] = FEMA_wrapper(fstem_imaging, fname_design, dirname_out, dirname_tabulated, dirname_imaging, datatype,...
'ranknorm', ranknorm, 'contrasts', contrasts, 'RandomEffects', RandomEffects, 'pihat_file', fname_pihat, 'nperms', nperms, 'mediation',mediation,'PermType',PermType,'tfce',tfce);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. Longitudinal with random effects F, A, S, E

fstem_imaging='longitudinal_FASE';
dirname_imaging = '/space/syn50/1/data/ABCD/d9smith/random_effects/behavioral/nih_tbx_longitudinal.txt';
% dirname_out = outDir; %filepath to save FEMA output

out_path = strcat(outDir, '/', designmat_file);
dirname_out = out_path;

RandomEffects = {'F','A','S','E'}; % Random effects to include: family, subject, error

% RUN FEMA
[fnames_out beta_hat beta_se zmat logpmat sig2tvec sig2mat beta_hat_perm beta_se_perm zmat_perm sig2tvec_perm sig2mat_perm logLikvec_perm inputs] = FEMA_wrapper(fstem_imaging, fname_design, dirname_out, dirname_tabulated, dirname_imaging, datatype,...
'ranknorm', ranknorm, 'contrasts', contrasts, 'RandomEffects', RandomEffects, 'pihat_file', fname_pihat, 'nperms', nperms, 'mediation',mediation,'PermType',PermType,'tfce',tfce);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5. Longitudinal with random effects F, A, D, S, E

fstem_imaging='longitudinal_FADSE';
dirname_imaging = '/space/syn50/1/data/ABCD/d9smith/random_effects/behavioral/nih_tbx_longitudinal.txt';
% dirname_out = outDir; %filepath to save FEMA output

out_path = strcat(outDir, '/', designmat_file);
dirname_out = out_path;

RandomEffects = {'F','A','D','S','E'}; % Random effects to include: family, subject, error

% RUN FEMA
[fnames_out beta_hat beta_se zmat logpmat sig2tvec sig2mat beta_hat_perm beta_se_perm zmat_perm sig2tvec_perm sig2mat_perm logLikvec_perm inputs] = FEMA_wrapper(fstem_imaging, fname_design, dirname_out, dirname_tabulated, dirname_imaging, datatype,...
'ranknorm', ranknorm, 'contrasts', contrasts, 'RandomEffects', RandomEffects, 'pihat_file', fname_pihat, 'nperms', nperms, 'mediation',mediation,'PermType',PermType,'tfce',tfce);

