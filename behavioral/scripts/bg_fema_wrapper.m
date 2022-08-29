%% Run FEMA for random effects project
%% Diana Smith
%% April 2022

% ADD CMIG tools directory to MATLAB path:
addpath(genpath('/home/d9smith/github/cmig_tools_internal'));

% Specify data release
dataRelease = '4.0';

% Path to store results
dirname_out = '/space/syn50/1/data/ABCD/d9smith/random_effects/behavioral/results/results_20220829';

% paths to inputs
measured_grm_file = fullfile('/space/amdale/1/tmp/ABCD_cache/abcd-sync/4.0/genomics/ABCD_rel4.0_grm.mat');
twin_grm_file = fullfile('/space/syn50/1/data/ABCD/d9smith/random_effects/behavioral/data/twins_assigned_grm.mat');
pheno_dir = '/space/syn50/1/data/ABCD/d9smith/random_effects/behavioral/data/pheno';

% Inputs that will remain the same for all models
atlasVersion = 'ABCD2_cor10';
dirname_tabulated = fullfile('/space/amdale/1/tmp/ABCD_cache/abcd-sync/4.0/tabulated/released/'); %KNOWN ISSUE: breaks when using txt files following NDA release --> must use pre-release csv files
fname_pregnancyID = fullfile('/home/sabad/requests/pregnancy_ID_07182022.csv');
fname_addressID = fullfile('/home/sabad/requests/recent_addr_07182022.csv'); 
fname_design = '/space/syn50/1/data/ABCD/d9smith/random_effects/behavioral/designMat/designMat0_empty.txt'; 

contrasts=[]; % Contrasts relate to columns in design matrix e.g. [1 -1] will take the difference between cols 1 and 2 in your design matrix (X)
ranknorm = 0; % Rank normalizes dependent variables (Y) (default = 0)
nperms = 0; % Number of permutations - set to 0 for now
mediation = 0; % If wanting to use outputs for a mediation analysis set mediation=1 - ensures same resampling scheme used for each model in fname_design
PermType = 'wildbootstrap'; %Default resampling method is null wild-bootstrap - to run mediation analysis need to use non-null wild-bootstrap ('wildboostrap-nn')
tfce = 0; % Columns in design matrix to loop over to calculate TFCE - selecting columns of interest improves efficiency - is `isempty(tfce_cols)` FEMA_wrapper will NOT run TFCE
RandomEstType = 'ML'; % specify random effects estimator (default is MoM)
Hessflag=0;
logLikflag=1;
ciflag=1;
datatype='external'; % can use txt with columns of ANY data type (e.g. ROIs, behavior) - runs mass univaraite LME across every column

%% Model 1: FAE Model, twins only at baseline, genetic relatedness assumed (will also run in OpenMx)
fstem_imaging = 'model1';
RandomEffects = {'F';'A';'E'}; 
fname_pihat = twin_grm_file; 
dirname_imaging = strcat(pheno_dir,'/','baseline_twins_res_agesexsite.txt');

% run FEMA
[fnames_out beta_hat beta_se zmat logpmat sig2tvec sig2mat beta_hat_perm beta_se_perm zmat_perm sig2tvec_perm sig2mat_perm inputs mask tfce_perm colnames_interest save_params logLikvec Hessmat] = FEMA_wrapper(fstem_imaging, fname_design, dirname_out, dirname_tabulated, dirname_imaging, datatype,...
'ranknorm', ranknorm, 'contrasts', contrasts, 'RandomEffects', RandomEffects, 'pihat_file', fname_pihat, 'nperms', nperms, 'mediation',mediation,'PermType',PermType,'tfce',tfce,'preg_file',fname_pregnancyID,'address_file',fname_addressID,...
'Hessflag',Hessflag,'ciflag',ciflag,'logLikflag',logLikflag,'RandomEstType',RandomEstType);

disp(sig2mat);
disp(logLikvec);

%% Model 2: FAE Model, twins only at baseline, with GRM included.
fstem_imaging = 
RandomEffects = {'F';'A';'E'};
fname_pihat = measured_grm_file;
dirname_imaging = 

% run FEMA
[fnames_out beta_hat beta_se zmat logpmat sig2tvec sig2mat beta_hat_perm beta_se_perm zmat_perm sig2tvec_perm sig2mat_perm inputs mask tfce_perm colnames_interest save_params logLikvec Hessmat] = FEMA_wrapper(fstem_imaging, fname_design, dirname_out, dirname_tabulated, dirname_imaging, datatype,...
'ranknorm', ranknorm, 'contrasts', contrasts, 'RandomEffects', RandomEffects, 'pihat_file', fname_pihat, 'nperms', nperms, 'mediation',mediation,'PermType',PermType,'tfce',tfce,'preg_file',fname_pregnancyID,'address_file',fname_addressID,...
'Hessflag',Hessflag,'ciflag',ciflag,'logLikflag',logLikflag,'RandomEstType',RandomEstType);

%% Model 3: FAE Model, full sample at baseline, with GRM included within family.
fstem_imaging = 
RandomEffects = {'F';'A';'E'};
fname_pihat = measured_grm_file;
dirname_imaging = 

% run FEMA
[fnames_out beta_hat beta_se zmat logpmat sig2tvec sig2mat beta_hat_perm beta_se_perm zmat_perm sig2tvec_perm sig2mat_perm inputs mask tfce_perm colnames_interest save_params logLikvec Hessmat] = FEMA_wrapper(fstem_imaging, fname_design, dirname_out, dirname_tabulated, dirname_imaging, datatype,...
'ranknorm', ranknorm, 'contrasts', contrasts, 'RandomEffects', RandomEffects, 'pihat_file', fname_pihat, 'nperms', nperms, 'mediation',mediation,'PermType',PermType,'tfce',tfce,'preg_file',fname_pregnancyID,'address_file',fname_addressID,...
'Hessflag',Hessflag,'ciflag',ciflag,'logLikflag',logLikflag,'RandomEstType',RandomEstType);

%% Model 4: FATE Model, full sample at baseline, with GRM within family.
fstem_imaging = 
RandomEffects = {'F';'A';'T';'E'};
fname_pihat = measured_grm_file;
dirname_imaging = 

% run FEMA
[fnames_out beta_hat beta_se zmat logpmat sig2tvec sig2mat beta_hat_perm beta_se_perm zmat_perm sig2tvec_perm sig2mat_perm inputs mask tfce_perm colnames_interest save_params logLikvec Hessmat] = FEMA_wrapper(fstem_imaging, fname_design, dirname_out, dirname_tabulated, dirname_imaging, datatype,...
'ranknorm', ranknorm, 'contrasts', contrasts, 'RandomEffects', RandomEffects, 'pihat_file', fname_pihat, 'nperms', nperms, 'mediation',mediation,'PermType',PermType,'tfce',tfce,'preg_file',fname_pregnancyID,'address_file',fname_addressID,...
'Hessflag',Hessflag,'ciflag',ciflag,'logLikflag',logLikflag,'RandomEstType',RandomEstType);

%% Model 5: FASTE Model, full sample at baseline and year 2, with GRM within family. 
% Note that all longitudinal analyses should include data that is preresidualized for age, sex, site, and practice effect.
fstem_imaging = 
RandomEffects = {'F';'A';'S';'T';'E'};
fname_pihat = measured_grm_file;
dirname_imaging = 

% run FEMA
[fnames_out beta_hat beta_se zmat logpmat sig2tvec sig2mat beta_hat_perm beta_se_perm zmat_perm sig2tvec_perm sig2mat_perm inputs mask tfce_perm colnames_interest save_params logLikvec Hessmat] = FEMA_wrapper(fstem_imaging, fname_design, dirname_out, dirname_tabulated, dirname_imaging, datatype,...
'ranknorm', ranknorm, 'contrasts', contrasts, 'RandomEffects', RandomEffects, 'pihat_file', fname_pihat, 'nperms', nperms, 'mediation',mediation,'PermType',PermType,'tfce',tfce,'preg_file',fname_pregnancyID,'address_file',fname_addressID,...
'Hessflag',Hessflag,'ciflag',ciflag,'logLikflag',logLikflag,'RandomEstType',RandomEstType);

%% Side question #1: discretizing zygosity
% Run models 3-5 with assigned zygosity.
fstem_imaging = 
RandomEffects = 
fname_pihat = 
dirname_imaging = 

% run FEMA
[fnames_out beta_hat beta_se zmat logpmat sig2tvec sig2mat beta_hat_perm beta_se_perm zmat_perm sig2tvec_perm sig2mat_perm inputs mask tfce_perm colnames_interest save_params logLikvec Hessmat] = FEMA_wrapper(fstem_imaging, fname_design, dirname_out, dirname_tabulated, dirname_imaging, datatype,...
'ranknorm', ranknorm, 'contrasts', contrasts, 'RandomEffects', RandomEffects, 'pihat_file', fname_pihat, 'nperms', nperms, 'mediation',mediation,'PermType',PermType,'tfce',tfce,'preg_file',fname_pregnancyID,'address_file',fname_addressID,...
'Hessflag',Hessflag,'ciflag',ciflag,'logLikflag',logLikflag,'RandomEstType',RandomEstType);

%% Side question #2: including fixed effect covariates
% Run models 1-5, preresidualized for all fixed effect covariates (genetic PCs, parental education, income).
% Note that only Model 1 will be reported in main text.
fstem_imaging = 
RandomEffects = 
fname_pihat = measured_grm_file;
dirname_imaging = 

% run FEMA
[fnames_out beta_hat beta_se zmat logpmat sig2tvec sig2mat beta_hat_perm beta_se_perm zmat_perm sig2tvec_perm sig2mat_perm inputs mask tfce_perm colnames_interest save_params logLikvec Hessmat] = FEMA_wrapper(fstem_imaging, fname_design, dirname_out, dirname_tabulated, dirname_imaging, datatype,...
'ranknorm', ranknorm, 'contrasts', contrasts, 'RandomEffects', RandomEffects, 'pihat_file', fname_pihat, 'nperms', nperms, 'mediation',mediation,'PermType',PermType,'tfce',tfce,'preg_file',fname_pregnancyID,'address_file',fname_addressID,...
'Hessflag',Hessflag,'ciflag',ciflag,'logLikflag',logLikflag,'RandomEstType',RandomEstType);

%% Side question #3: are twins necessary?
% Run FASE model, full sample at baseline and year 2, with GRM within family.
% Note that all longitudinal analyses should include data that is preresidualized for age, sex, site, and practice effect.
fstem_imaging = 
RandomEffects = 
fname_pihat = measured_grm_file; 
dirname_imaging =

% run FEMA
[fnames_out beta_hat beta_se zmat logpmat sig2tvec sig2mat beta_hat_perm beta_se_perm zmat_perm sig2tvec_perm sig2mat_perm inputs mask tfce_perm colnames_interest save_params logLikvec Hessmat] = FEMA_wrapper(fstem_imaging, fname_design, dirname_out, dirname_tabulated, dirname_imaging, datatype,...
'ranknorm', ranknorm, 'contrasts', contrasts, 'RandomEffects', RandomEffects, 'pihat_file', fname_pihat, 'nperms', nperms, 'mediation',mediation,'PermType',PermType,'tfce',tfce,'preg_file',fname_pregnancyID,'address_file',fname_addressID,...
'Hessflag',Hessflag,'ciflag',ciflag,'logLikflag',logLikflag,'RandomEstType',RandomEstType);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% extra code, no longer used
if 0
    % paste extra code here
    % specify array of random effects
    random_effects = {{'F','A','E'};{'F','A','T','E'};{'F','A','T','H','E'};{'F','A','T','S','E'};{'F','A','T','H','S','E'};
    {'H','A','E'};{'H','A','T','E'};{'H','A','T','S','E'}};

    % specify array of design 
    % fname_design = '/home/d9smith/projects/random_effects/behavioral/designMat/designMat1_allcovs.txt'; % for debugging only
    designmat_dir = '/home/d9smith/projects/random_effects/behavioral/designMat';
    designmat_file = dir(sprintf('%s/designMat*.txt', designmat_dir));
    designmat_file = {designmat_file.name}';
    designmat_array = strcat(designmat_dir, '/', designmat_file);

    % when only running one model
    %fname_design = '/space/syn50/1/data/ABCD/d9smith/random_effects/behavioral/designMat/designMat0_empty.txt' % baseline only
    fname_design = '/space/syn50/1/data/ABCD/d9smith/random_effects/behavioral/designMat/designMat1_allcovs.txt'
    RandomEffects = {'F';'A';'S';'T';'E'};
    dirname_imaging = '/space/syn50/1/data/ABCD/d9smith/random_effects/behavioral/data/pheno/baseline_twins_res_agesexsite.txt'; 
    fstem_imaging = sprintf('%s',RandomEffects{:});
    dirname_out = strcat(outDir, '/',dirname_imaging(71:end-4)); % when looping through design matrices
    % dirname_out = strcat(outDir, '/',fname_design(1:-4)); % when only using one design matrix
    disp(dirname_out);

    % loop through each design matrix
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
            [fnames_out beta_hat beta_se zmat logpmat sig2tvec sig2mat beta_hat_perm beta_se_perm zmat_perm sig2tvec_perm sig2mat_perm inputs mask tfce_perm colnames_interest save_params logLikvec Hessmat] = FEMA_wrapper(fstem_imaging, fname_design, dirname_out, dirname_tabulated, dirname_imaging, datatype,...
            'ranknorm', ranknorm, 'contrasts', contrasts, 'RandomEffects', RandomEffects, 'pihat_file', fname_pihat, 'nperms', nperms, 'mediation',mediation,'PermType',PermType,'tfce',tfce,'preg_file',fname_pregnancyID,'address_file',fname_addressID,...
            'Hessflag',Hessflag,'ciflag',ciflag,'logLikflag',logLikflag,'RandomEstType',RandomEstType);
            
            % save sig2mat_perm as output
            save(fnames_out{:}, 'sig2mat_perm', '-append');
            save(fnames_out{:}, 'logLikvec_perm', '-append');
            %writematrix(sig2mat_perm, sprintf('%s_sig2mat_perm.csv',fnames_out{:}(1:end-4)))
            writematrix(sig2mat,sprintf('%s/sig2mat_%s.csv',dirname_out,[random_effects{j}{:}]));
        end
    end
end