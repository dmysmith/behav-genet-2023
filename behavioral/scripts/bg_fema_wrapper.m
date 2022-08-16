%% Run FEMA for random effects project
%% Diana Smith
%% April 2022

% - add simulated phenotype

% Note: all models adjust for age and sex.

% Model 1: (will run in OpenMx and FEMA)
% FAE Model, twins only at baseline, genetic relatedness assumed
% Model 2: (will run in OpenMx and FEMA)
% FAE Model, twins only at baseline, with GRM included.
% Model 3: 
% 3a: FAE Model, full sample at baseline, with GRM included for everyone.
% 3b: FAE Model, full sample at baseline, with GRM included only within family.
% Model 4:
% FAE Model, full sample at baseline, with GRM within family, plus covariates (genetic PCs, parental education, income).
% Model 5:
% 5a: FATE Model, full sample at baseline, with GRM (will explain more variance than prev model) +/- fixed covariates
% 5b: FATHE Model, full sample at baseline, with GRM (will explain more variance than prev model) +/- fixed covariates 
% Model 6: 
% FATSE Model, full sample at two timepoints, with GRM +/- covariates (maybe add address effect to this one?)

% ADD CMIG tools directory to MATLAB path:
addpath(genpath('/home/d9smith/github/cmig_tools_internal'));

% Specify data release
dataRelease = '4.0';

% Path to store results
outDir = '/space/syn50/1/data/ABCD/d9smith/random_effects/behavioral';

% specify array of random effects
random_effects = {{'F','A','E'};{'F','A','T','E'};{'F','A','T','H','E'};{'F','A','T','S','E'};{'F','A','T','H','S','E'}};

% specify array of design 
% fname_design = '/home/d9smith/projects/random_effects/behavioral/designMat/designMat1_allcovs.txt'; % for debugging only
designmat_dir = '/home/d9smith/projects/random_effects/behavioral/designMat';
designmat_file = dir(sprintf('%s/designMat*.txt', designmat_dir));
designmat_file = {designmat_file.name}';
designmat_array = strcat(designmat_dir, '/', designmat_file);

% specify imaging file (the outcome you are predicting)
dirname_imaging_baseline = '/space/syn50/1/data/ABCD/d9smith/random_effects/behavioral/nih_tbx_baseline.txt';
dirname_imaging_longitudinal = '/space/syn50/1/data/ABCD/d9smith/random_effects/behavioral/nih_tbx_longitudinal.txt';

% Inputs to FEMA_wrapper.m
atlasVersion = 'ABCD2_cor10';
dirname_tabulated = fullfile('/space/amdale/1/tmp/ABCD_cache/abcd-sync/4.0/tabulated/released/'); %KNOWN ISSUE: breaks when using txt files following NDA release --> must use pre-release csv files
fname_pihat = fullfile('/space/amdale/1/tmp/ABCD_cache/abcd-sync/4.0/genomics/ABCD_rel4.0_grm.mat');
fname_addressID = fullfile('/home/sabad/requests/recent_addr_07182022.csv'); 
fname_pregnancyID = fullfile('/home/sabad/requests/pregnancy_ID_07182022.csv');

% Optional inputs for `FEMA_wrapper.m` depending on analysis
contrasts=[]; % Contrasts relate to columns in design matrix e.g. [1 -1] will take the difference between cols 1 and 2 in your design matrix (X)
ranknorm = 0; % Rank normalizes dependent variables (Y) (default = 0)
nperms = 0; % Number of permutations - set to 0 for now
mediation = 0; % If wanting to use outputs for a mediation analysis set mediation=1 - ensures same resampling scheme used for each model in fname_design
PermType = 'wildbootstrap'; %Default resampling method is null wild-bootstrap - to run mediation analysis need to use non-null wild-bootstrap ('wildboostrap-nn')
tfce = 0; % Columns in design matrix to loop over to calculate TFCE - selecting columns of interest improves efficiency - is `isempty(tfce_cols)` FEMA_wrapper will NOT run TFCE

datatype='external'; % can use txt with columns of ANY data type (e.g. ROIs, behavior) - runs mass univaraite LME across every column

%% loop through each design matrix
for i = 1:numel(designmat_array)
   fname_design = designmat_array(i);
   %% loop through each set of random effects
   for j = 1:numel(random_effects)     
        % define dynamic inputs
        fstem_imaging = sprintf('%s',random_effects{j}{:});
        dirname_out = strcat(outDir, '/',designmat_file{i}(1:end-4)); % when looping through design matrices
        % dirname_out = strcat(outDir, '/',fname_design(1:-4)); % when only using one design matrix
        disp(dirname_out);
        RandomEffects = random_effects{j};
        disp(RandomEffects);

        % for models including S, we want to use the longitudinal outcome data, otherwise only baseline
        if any(strcmp(random_effects{j},'S'))
            dirname_imaging = dirname_imaging_longitudinal;
        else
            dirname_imaging = dirname_imaging_baseline;
        end
        
        % run FEMA
        [fnames_out beta_hat beta_se zmat logpmat sig2tvec sig2mat beta_hat_perm beta_se_perm zmat_perm sig2tvec_perm sig2mat_perm logLikvec_perm inputs] = FEMA_wrapper(fstem_imaging, fname_design, dirname_out, dirname_tabulated, dirname_imaging, datatype,...
        'ranknorm', ranknorm, 'contrasts', contrasts, 'RandomEffects', RandomEffects, 'pihat_file', fname_pihat, 'nperms', nperms, 'mediation',mediation,'PermType',PermType,'tfce',tfce,'preg_file',fname_pregnancyID,'address_file',fname_addressID);
        
        % save sig2mat_perm as output
        save(fnames_out{:}, 'sig2mat_perm', '-append');
        save(fnames_out{:}, 'logLikvec_perm', '-append');
        %writematrix(sig2mat_perm, sprintf('%s_sig2mat_perm.csv',fnames_out{:}(1:end-4)))
        writematrix(sig2mat,sprintf('%s/sig2mat_%s.csv',dirname_out,[random_effects{j}{:}]));
    end
end

