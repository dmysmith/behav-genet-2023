%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Export FEMA random effects results to csv
%% Diana Smith
%% March 2022

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


load('/space/syn50/1/data/ABCD/d9smith/reading/agesex_only/FEMA_results_external_baseline_FE.mat')  
writematrix(sig2mat,'/home/d9smith/projects/reading/agesex_only/sig2mat_baseline_FE.csv')  

load('/space/syn50/1/data/ABCD/d9smith/reading/agesex_only/FEMA_results_external_baseline_FAE.mat')
writematrix(sig2mat,'/home/d9smith/projects/reading/agesex_only/sig2mat_baseline_FAE.csv')  

load('/space/syn50/1/data/ABCD/d9smith/reading/agesex_only/FEMA_results_external_longitudinal_FSE.mat')
writematrix(sig2mat,'/home/d9smith/projects/reading/agesex_only/sig2mat_longitudinal_FSE.csv')  

load('/space/syn50/1/data/ABCD/d9smith/reading/agesex_only/FEMA_results_external_longitudinal_FASE.mat')
writematrix(sig2mat,'/home/d9smith/projects/reading/agesex_only/sig2mat_longitudinal_FASE.csv')  