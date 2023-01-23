%% Run FEMA for behavior genetics paper 
%% Diana Smith
%% April 2022

% ADD CMIG tools directory to MATLAB path:
% Note that I used the internal repo (not the public one) 
% because we added an option to use a maximum likelihood estimator.
addpath(genpath('/home/d9smith/github/cmig_tools_internal'));
% addpath(genpath('/home/dale/matlab/FEMA'));

% Specify data release
dataRelease = '4.0';

% Path to store results
dirname_out = '/space/syn50/1/data/ABCD/d9smith/random_effects/behavioral/results/results_20230107';

% paths to inputs
measured_grm_file = fullfile('/space/amdale/1/tmp/ABCD_cache/abcd-sync/4.0/genomics/ABCD_rel4.0_grm.mat');
assigned_grm_file = fullfile('/space/syn50/1/data/ABCD/d9smith/random_effects/behavioral/data/all_discrete_grm.mat');
twin_grm_file = fullfile('/space/syn50/1/data/ABCD/d9smith/random_effects/behavioral/data/twins_assigned_grm.mat');
pheno_dir = '/space/syn50/1/data/ABCD/d9smith/random_effects/behavioral/data/pheno';

% Inputs that will remain the same for all models
atlasVersion = 'ABCD2_cor10';
dirname_tabulated = fullfile('/space/amdale/1/tmp/ABCD_cache/abcd-sync/4.0/tabulated/released/'); %KNOWN ISSUE: breaks when using txt files following NDA release --> must use pre-release csv files
fname_pregnancyID = fullfile('/home/sabad/requests/pregnancy_ID_01172023.csv');
fname_addressID = fullfile('/home/sabad/requests/recent_addr_07182022.csv'); 
fname_design = '/space/syn50/1/data/ABCD/d9smith/random_effects/behavioral/designMat/designMat0_empty.txt'; 

contrasts=[]; % Contrasts relate to columns in design matrix e.g. [1 -1] will take the difference between cols 1 and 2 in your design matrix (X)
ranknorm = 0; % Rank normalizes dependent variables (Y) (default = 0)
nperms = 0; % Number of permutations - set to 0 for now
mediation = 0; % If wanting to use outputs for a mediation analysis set mediation=1 - ensures same resampling scheme used for each model in fname_design
PermType = 'wildbootstrap-nn'; %Default resampling method is null wild-bootstrap - to run mediation analysis need to use non-null wild-bootstrap ('wildboostrap-nn')
tfce = 0; % Columns in design matrix to loop over to calculate TFCE - selecting columns of interest improves efficiency - is `isempty(tfce_cols)` FEMA_wrapper will NOT run TFCE
RandomEstType = 'ML'; % specify random effects estimator (default is MoM)
Hessflag=0;
logLikflag=1;
ciflag=1;
datatype='external'; % can use txt with columns of ANY data type (e.g. ROIs, behavior) - runs mass univaraite LME across every column

fstem_imaging = {}; RandomEffects = {}; fname_pihat = {}; dirname_imaging = {}; titles = {};

%% Model 1: ACE Model, twins only at baseline, genetic relatedness assumed (will also run in OpenMx)
i=1;
fstem_imaging{i} = 'model1';
titles{i} = 'FEMA ACE Model, twin sample, baseline, kinship-derived GRM';
RandomEffects{i} = {'A';'F';'E'}; 
fname_pihat{i} = twin_grm_file; 
dirname_imaging{i} = strcat(pheno_dir,'/','baseline_twins_res_agesex.txt');

%% Model 2: ACE Model, twins only at baseline, with GRM included.
i=2;
fstem_imaging{i} = 'model2';
titles{i} =  'ACE Model, twin sample, baseline, SNP-derived GRM';
RandomEffects{i} = {'A';'F';'E'};
fname_pihat{i} = measured_grm_file;
dirname_imaging{i} = strcat(pheno_dir,'/','baseline_twins_res_agesex.txt'); 

%% Model 3: ACE Model, full sample at baseline, with GRM included within family.
i=3;
fstem_imaging{i} = 'model3';
titles{i} = 'ACE Model, full sample, baseline, SNP-derived GRM';
RandomEffects{i} = {'A';'F';'E'};
fname_pihat{i} = measured_grm_file;
dirname_imaging{i} = strcat(pheno_dir,'/','baseline_full_res_agesex.txt');  

%% Model 4: ACTE Model, full sample at baseline, with GRM within family.
i=4;
fstem_imaging{i} = 'model4';
titles{i} = 'ACTE Model, full sample, baseline, SNP-derived GRM';
RandomEffects{i} = {'A';'F';'T';'E'};
fname_pihat{i} = measured_grm_file;
dirname_imaging{i} = strcat(pheno_dir,'/','baseline_full_res_agesex.txt'); 

%% Model 5: ACTE Model, full sample at baseline and year 2, with GRM within family. 
% Note that all longitudinal analyses should include data that is preresidualized for age, sex, site, and practice effect.
i=5;
fstem_imaging{i} = 'model5';
titles{i} = 'ACTE Model, full sample, baseline and Y2, SNP-derived GRM';
RandomEffects{i} = {'A';'F';'T';'E'};
fname_pihat{i} = measured_grm_file;
dirname_imaging{i} = strcat(pheno_dir,'/','longitudinal_full_res_agesexprac.txt');   

%% Model 6: ACTSE Model, full sample at baseline and year 2, with GRM within family. 
% Note that all longitudinal analyses should include data that is preresidualized for age, sex, site, and practice effect.
i=6;
fstem_imaging{i} = 'model6';
titles{i} = 'ACTSE Model, full sample, baseline and Y2, SNP-derived GRM';
RandomEffects{i} = {'A';'F';'T';'S';'E'};
fname_pihat{i} = measured_grm_file;
dirname_imaging{i} = strcat(pheno_dir,'/','longitudinal_full_res_agesexprac.txt');   

%% Side question #1: discretizing zygosity
% Run models 3-5 with assigned zygosity.
i=7;
fstem_imaging{i} = 's1_assigngrm_m3';
titles{i} = 'ACE Model, full sample, baseline, kinship-derived GRM';
RandomEffects{i} = {'A';'F';'E'};
fname_pihat{i} = assigned_grm_file;
dirname_imaging{i} = strcat(pheno_dir,'/','baseline_full_res_agesex.txt');  

i=8;
fstem_imaging{i} = 's1_assigngrm_m4';
titles{i} = 'ACTE Model, full sample, baseline, kinship-derived GRM';
RandomEffects{i} = {'A';'F';'T';'E'}; 
fname_pihat{i} = assigned_grm_file; 
dirname_imaging{i} = strcat(pheno_dir,'/','baseline_full_res_agesex.txt');  

i=9;
fstem_imaging{i} = 's1_assigngrm_m5';
titles{i} = 'ACTSE Model, full sample, baseline and Y2, kinship-derived GRM';
RandomEffects{i} = {'A';'F';'T';'S';'E'};
fname_pihat{i} = assigned_grm_file;
dirname_imaging{i} = strcat(pheno_dir,'/','longitudinal_full_res_agesexprac.txt'); 

%% Side question #2: including fixed effect covariates
% Run models 1-5, preresidualized for all fixed effect covariates (genetic PCs, parental education, income).
% Note that only Model 1 will be reported in main text.
i=10;
fstem_imaging{i} = 's2_allcovs_m1';
titles{i} = 'ACE Model, twin sample, baseline, kinship-derived GRM, all covariates';
RandomEffects{i} = {'A';'F';'E'}; 
fname_pihat{i} = twin_grm_file;
dirname_imaging{i} = strcat(pheno_dir,'/','baseline_twins_res_agesexsiteeducincpcs.txt'); 

i=11;
fstem_imaging{i} = 's2_allcovs_m2';
titles{i} = 'ACE Model, twin sample, baseline, SNP-derived GRM, all covariates';
RandomEffects{i} = {'A';'F';'E'}; 
fname_pihat{i} = measured_grm_file;
dirname_imaging{i} = strcat(pheno_dir,'/','baseline_twins_res_agesexsiteeducincpcs.txt');

i=12;
fstem_imaging{i} = 's2_allcovs_m3';
titles{i} = 'ACE Model, full sample, baseline, SNP-derived GRM, all covariates';
RandomEffects{i} = {'A';'F';'E'};
fname_pihat{i} = measured_grm_file;
dirname_imaging{i} = strcat(pheno_dir,'/','baseline_full_res_agesexsiteeducincpcs.txt');

i=13;
fstem_imaging{i} = 's2_allcovs_m4';
titles{i} = 'ACTE Model, full sample, baseline, SNP-derived GRM, all covariates';
RandomEffects{i} = {'A';'F';'T';'E'}; 
fname_pihat{i} = measured_grm_file;
dirname_imaging{i} = strcat(pheno_dir,'/','baseline_twins_res_agesexsiteeducincpcs.txt');

i=14;
fstem_imaging{i} = 's2_allcovs_m5';
titles{i} = 'ACTSE Model, full sample, baseline and Y2, SNP-derived GRM, all covariates';
RandomEffects{i} = {'A';'F';'T';'S';'E'}; 
fname_pihat{i} = measured_grm_file;
dirname_imaging{i} = strcat(pheno_dir,'/','longitudinal_full_res_agesexsitepraceducincpcs.txt');

%% Side question #3: are twins necessary?
% Run FASE (ACSE) model, full sample minus twins at baseline and year 2, with GRM within family.
% Note that all longitudinal analyses should include data that is preresidualized for age, sex, site, and practice effect.
i=15;
fstem_imaging{i} = 's3_notwins_comparison';
titles{i} = 'ACSE model, full sample, baseline and Y2, SNP-derived GRM';
RandomEffects{i} = {'A';'F';'S';'E';}
fname_pihat{i} = measured_grm_file; 
dirname_imaging{i} = strcat(pheno_dir,'/','longitudinal_full_res_agesexprac.txt');

i=16;
fstem_imaging{i} = 's3_notwins';
titles{i} = 'ACSE model, full sample minus twin sample, baseline and Y2, SNP-derived GRM';
RandomEffects{i} = {'A';'F';'S';'E';}
fname_pihat{i} = measured_grm_file; 
dirname_imaging{i} = strcat(pheno_dir,'/','longitudinal_notwins_res_agesexprac.txt');

if 0
  %% Side Question #4: What is the deal with confidence intervals?
  nperms = 1000;

  i=17;
  fstem_imaging{i} = 's4_1000perms_model1';
  titles{i} = 'FEMA ACE Model, twin sample, baseline, assigned GRM, 1000 permutations';
  RandomEffects{i} = {'A';'F';'E'}; 
  fname_pihat{i} = twin_grm_file; 
  dirname_imaging{i} = strcat(pheno_dir,'/','baseline_twins_res_agesex.txt');

  i=18;
  fstem_imaging{i} = 's4_1000perms_model2';
  titles{i} =  'ACE Model, twin sample, baseline, measured GRM, 1000 permutations';
  RandomEffects{i} = {'A';'F';'E'};
  fname_pihat{i} = measured_grm_file;
  dirname_imaging{i} = strcat(pheno_dir,'/','baseline_twins_res_agesex.txt'); 

  i=19;
  fstem_imaging{i} = 's4_1000perms_model3';
  titles{i} = 'ACE Model, full sample, baseline, measured GRM, 1000 permutations';
  RandomEffects{i} = {'A';'F';'E'};
  fname_pihat{i} = measured_grm_file;
  dirname_imaging{i} = strcat(pheno_dir,'/','baseline_full_res_agesex.txt');  

  i=20;
  fstem_imaging{i} = 's4_1000perms_model4';
  titles{i} = 'ACTE Model, full sample, baseline, measured GRM, 1000 permutations';
  RandomEffects{i} = {'A';'F';'T';'E'};
  fname_pihat{i} = measured_grm_file;
  dirname_imaging{i} = strcat(pheno_dir,'/','baseline_full_res_agesex.txt'); 

  % Note that all longitudinal analyses should include data that is preresidualized for age, sex, site, and practice effect.
  i=21;
  fstem_imaging{i} = 's4_1000perms_model6';
  titles{i} = 'ACTSE Model, full sample, baseline and Y2, measured GRM, 1000 permutations';
  RandomEffects{i} = {'A';'F';'T';'S';'E'};
  fname_pihat{i} = measured_grm_file;
  dirname_imaging{i} = strcat(pheno_dir,'/','longitudinal_full_res_agesexprac.txt');   

  % i=17;
  % fstem_imaging{i} = 'test_y2';
  % titles{i} = 'ACE model, full sample at year 2 only, with GRM, age, sex';
  % RandomEffects{i} = {'A';'F';'E'};
  % fname_pihat{i} = measured_grm_file; 
  % dirname_imaging{i} = strcat(pheno_dir,'/','y2_full_res_agesex.txt');

  % i=18;
  % fstem_imaging{i} = 'test_y2_allcovs';
  % titles{i} = 'ACE model, full sample at year 2 only, measured GRM, all covariates';
  % RandomEffects{i} = {'A';'F';'E'};
  % fname_pihat{i} = measured_grm_file; 
  % dirname_imaging{i} = strcat(pheno_dir,'/','y2_full_res_agesexsiteeducincpcs.txt')d;
end

save(strcat(dirname_out, '/', 'model_parameters.mat'), 'fstem_imaging', 'titles', 'RandomEffects');

% run FEMA
for i = 1:length(fstem_imaging)
% for i = [16,17]
    if isempty(fstem_imaging{i})
        continue;
    end
    
    [fnames_out beta_hat beta_se zmat logpmat sig2tvec sig2mat beta_hat_perm beta_se_perm zmat_perm sig2tvec_perm sig2mat_perm inputs mask tfce_perm colnames_interest save_params logLikvec Hessmat] = FEMA_wrapper(fstem_imaging{i}, fname_design, dirname_out, dirname_tabulated, dirname_imaging{i}, datatype,...
    'ranknorm', ranknorm, 'contrasts', contrasts, 'RandomEffects', RandomEffects{i}, 'pihat_file', fname_pihat{i}, 'nperms', nperms, 'mediation',mediation,'PermType',PermType,'tfce',tfce,'preg_file',fname_pregnancyID,'address_file',fname_addressID,...
    'Hessflag',Hessflag,'ciflag',ciflag,'logLikflag',logLikflag,'RandomEstType',RandomEstType);
    save(fnames_out{:}, 'logLikvec', '-append');
    save(fnames_out{:}, 'fstem_imaging', '-append');
    save(fnames_out{:}, 'RandomEffects', '-append');

    json_out = sprintf('%s.json',fnames_out{:}(1:end-4));
    sig2mat_out = sprintf('%s_sig2mat.json',fnames_out{:}(1:end-4));
    
    fid = fopen(json_out,'w');
    s = struct('beta_hat', beta_hat, 'beta_se', beta_se, 'zmat', zmat, 'logpmat', logpmat, 'sig2tvec', sig2tvec, ...
      'colnames_interest', colnames_interest, 'logLikvec', logLikvec);
    fprintf(fid,jsonencode(s));
    fclose(fid);

    fid = fopen(sig2mat_out,'w');
    s = struct('sig2mat', sig2mat); 
    % fprintf(fid,jsonencode(s));
    fprintf(fid,jsonencode(sig2mat));  
    fclose(fid);

end