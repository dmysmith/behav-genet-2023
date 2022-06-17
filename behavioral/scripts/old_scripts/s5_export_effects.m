%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Export FEMA random effects results to csv
%% Diana Smith
%% April 2022

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load('/space/syn50/1/data/ABCD/d9smith/random_effects/behavioral/FEMA_wrapper_output_external_baseline_FE.mat')  
% writematrix(sig2mat,'/home/d9smith/projects/random_effects/behavioral/sig2mat_baseline_FE.csv')  
if 0
    load('/space/syn50/1/data/ABCD/d9smith/random_effects/behavioral/FEMA_wrapper_output_external_baseline_FAE.mat');
    writematrix(sig2mat,'/home/d9smith/projects/random_effects/behavioral/sig2mat_baseline_FAE.csv');
    
    load('/space/syn50/1/data/ABCD/d9smith/random_effects/behavioral/FEMA_wrapper_output_external_longitudinal_FSE.mat');
    writematrix(sig2mat,'/home/d9smith/projects/random_effects/behavioral/sig2mat_longitudinal_FSE.csv');
    
    load('/space/syn50/1/data/ABCD/d9smith/random_effects/behavioral/FEMA_wrapper_output_external_longitudinal_FASE.mat');
    writematrix(sig2mat,'/home/d9smith/projects/random_effects/behavioral/sig2mat_longitudinal_FASE.csv');
end

load('/space/syn50/1/data/ABCD/d9smith/random_effects/behavioral/designMat1_allcovs/nullWB_1000perms/FEMA_wrapper_output_external_FADSE.mat');
writematrix(sig2mat_perm,'/space/syn50/1/data/ABCD/d9smith/random_effects/behavioral/designMat1_allcovs/nullWB_1000perms/FEMA_wrapper_output_external_FADSE_sig2mat_perm.csv');

load('/space/syn50/1/data/ABCD/d9smith/random_effects/behavioral/designMat1_allcovs/nullWB_1000perms/FEMA_wrapper_output_external_FASE.mat');
writematrix(sig2mat_perm,'/space/syn50/1/data/ABCD/d9smith/random_effects/behavioral/designMat1_allcovs/nullWB_1000perms/FEMA_wrapper_output_external_FASE_sig2mat_perm.csv');

load('/space/syn50/1/data/ABCD/d9smith/random_effects/behavioral/designMat1_allcovs/nullWB_1000perms/FEMA_wrapper_output_external_FSE.mat');
writematrix(sig2mat_perm,'/space/syn50/1/data/ABCD/d9smith/random_effects/behavioral/designMat1_allcovs/nullWB_1000perms/FEMA_wrapper_output_external_FSE_sig2mat_perm.csv');

load('/space/syn50/1/data/ABCD/d9smith/random_effects/behavioral/designMat1_allcovs/nullWB_1000perms/FEMA_wrapper_output_external_FAE.mat');
writematrix(sig2mat_perm,'/space/syn50/1/data/ABCD/d9smith/random_effects/behavioral/designMat1_allcovs/nullWB_1000perms/FEMA_wrapper_output_external_FAE_sig2mat_perm.csv');