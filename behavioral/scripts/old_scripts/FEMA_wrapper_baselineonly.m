%% Run FEMA for random effects project
%% Diana Smith
%% April 2022

% ADD CMIG tools directory to MATLAB path:
addpath(genpath('/home/d9smith/github/cmig_tools_internal'));

% Specify data release
dataRelease = '4.0';

% Path to store results
outDir = '/space/syn50/1/data/ABCD/d9smith/random_effects/behavioral/baseline';

% specify array of random effects
random_effects = {{'F','A','E'}};

% specify array of design 
% fname_design = '/home/d9smith/projects/random_effects/behavioral/designMat1_allcovs.txt'; 
designmat_dir = '/home/d9smith/projects/random_effects/behavioral';
designmat_file = dir(sprintf('%s/designMat*.txt', designmat_dir));
designmat_file = {designmat_file.name}';
designmat_array = strcat(designmat_dir, '/', designmat_file);

% specify imaging file (the outcome you are predicting)
dirname_imaging = '/space/syn50/1/data/ABCD/d9smith/random_effects/behavioral/nih_tbx_baseline.txt';

% Inputs to FEMA_wrapper.m
atlasVersion = 'ABCD2_cor10';
dirname_tabulated = fullfile('/space/amdale/1/tmp/ABCD_cache/abcd-sync/4.0/tabulated/released/'); %KNOWN ISSUE: breaks when using txt files following NDA release --> must use pre-release csv files
% fname_pihat = fullfile(abcd_sync_path, '3.0','genomics','ABCD_rel3.0_pihat.mat'); %FIXME: replace with 4.0 when available
fname_pihat = fullfile('/space/amdale/1/tmp/ABCD_cache/abcd-sync/3.0/genomics/ABCD_rel3.0_pihat.mat');

% Optional inputs for `FEMA_wrapper.m` depending on analysis
contrasts=[]; % Contrasts relate to columns in design matrix e.g. [1 -1] will take the difference between cols 1 and 2 in your design matrix (X)
ranknorm = 0; % Rank normalizes dependent variables (Y) (default = 0)
nperms = 10000; % Number of permutations - if wanting to use resampling methods nperms>0 and load Anders' script for FEMA_fit
mediation = 0; % If wanting to use outputs for a mediation analysis set mediation=1 - ensures same resampling scheme used for each model in fname_design
PermType = 'wildbootstrap'; %Default resampling method is null wild-bootstrap - to run mediation analysis need to use non-null wild-bootstrap ('wildboostrap-nn')
tfce = 0; % Columns in design matrix to loop over to calculate TFCE - selecting columns of interest improves efficiency - is `isempty(tfce_cols)` FEMA_wrapper will NOT run TFCE

datatype='external'; % can use txt with columns of ANY data type (e.g. ROIs, behavior) - runs mass univaraite LME across every column

% if nperms > 0 will need to add Anders' script
% if nperms > 0
if 0
    addpath(genpath('~dale/matlab/SSE_local'));
    % outDir = sprintf('%s/%.0fperms', outDir, nperms);
end

%% loop through each design matrix
for i = 1:numel(designmat_array)
   fname_design = designmat_array(i);
   %% loop through each set of random effects
   for j = 1:numel(random_effects)     
        % define dynamic inputs
        fstem_imaging = sprintf('%s',random_effects{j}{:});
        dirname_out = strcat(outDir, '/',designmat_file{i}(1:end-4)); % when looping through design matrices
        % dirname_out = strcat(outDir, '/',designmat_file(1:end-4)); % when only using one design matrix
        disp(dirname_out);
        RandomEffects = random_effects{j};
        disp(RandomEffects);
        
        % run FEMA
        [fnames_out beta_hat beta_se zmat logpmat sig2tvec sig2mat beta_hat_perm beta_se_perm zmat_perm sig2tvec_perm sig2mat_perm logLikvec_perm inputs] = FEMA_wrapper(fstem_imaging, fname_design, dirname_out, dirname_tabulated, dirname_imaging, datatype,...
        'ranknorm', ranknorm, 'contrasts', contrasts, 'RandomEffects', RandomEffects, 'pihat_file', fname_pihat, 'nperms', nperms, 'mediation',mediation,'PermType',PermType,'tfce',tfce);
        
        % save sig2mat_perm as output
        save(fnames_out{:}, 'sig2mat_perm', '-append');
        save(fnames_out{:}, 'logLikvec_perm', '-append');
        %writematrix(sig2mat_perm, sprintf('%s_sig2mat_perm.csv',fnames_out{:}(1:end-4)))
        writematrix(sig2mat,sprintf('%s/sig2mat_%s.csv',dirname_out,[random_effects{j}{:}]));
    end
end
% export effects
