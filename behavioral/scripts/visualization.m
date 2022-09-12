% visualization of variances
% diana smith
% may 2022

% inputs
results_dir = "/space/syn50/1/data/ABCD/d9smith/random_effects/behavioral/results/results_20220901";
outpath = "/home/d9smith/projects/random_effects/behavioral/results/plots"; % outpath = "/home/d9smith/tmp";

addpath(genpath('/home/d9smith/projects/random_effects/behavioral/scripts'));

param_file = strcat(results_dir, '/', 'model_parameters.mat');
load(param_file);

results_file = strcat(results_dir,"/FEMA_wrapper_output_external_",fstem_imaging',".mat");

% for i=1:size(fstem_imaging,2)
for i=[16,17]
    visualize_ds(results_file{i}, outpath, fstem_imaging, titles, i)
end

% load openmx data

mxa_file = strcat(results_dir,"/openmx_A.csv"); 
mxc_file = strcat(results_dir,"/openmx_C.csv"); 
mxe_file = strcat(results_dir,"/openmx_E.csv"); 
mxloglik_file = strcat(results_dir,"/openmx_loglik.csv"); 

mxa = readtable(mxa_file);
mxc = readtable(mxc_file);
mxe = readtable(mxe_file);
mxloglik = readtable(mxloglik_file);

i=15;
fstem_imaging{i} = "openmx";
titles{i} = "OpenMx ACE/FAE Model, twins only at baseline, GRM assumed";
RandomEffects{i} = {'F';'A';'E'};
colnames_imaging = mxa.task';

sig2mat(:,:,1) = [mxc.openmx';mxa.openmx';mxe.openmx'];
sig2mat(:,:,2) = [mxc.openmx_ci_lower';mxa.openmx_ci_lower';mxe.openmx_ci_lower'];
sig2mat(:,:,3) = [mxc.openmx_ci_upper';mxa.openmx_ci_upper';mxe.openmx_ci_upper'];

mx_results_file = strcat(results_dir, '/', 'openmx.mat');

save(mx_results_file, 'sig2mat','RandomEffects','colnames_imaging');

visualize_ds(mx_results_file, outpath, fstem_imaging, titles, i);

if 0
    fname_design = 'designMat1_allcovs';
    nperms = 1000;
    RandomEffects = "FASE";
    inpath = "/space/syn50/1/data/ABCD/d9smith/random_effects/behavioral";

    results_file = sprintf("%s/%s/nullWB_%dperms/FEMA_wrapper_output_external_%s.mat",inpath,fname_design,nperms,RandomEffects);

    % calculate percentiles
    perc = nan(size(sig2mat));
    for i = 1:size(sig2mat,1)
        for j = 1:size(sig2mat,2)
            perc(i,j,1) = prctile(sig2mat_perm(i,j,2:end), 2.5);
            perc(i,j,2) = prctile(sig2mat_perm(i,j,2:end), 97.5);
        end
    end
    
    title(sprintf('Model including: %s, %s',fname_design(12:end),RandomEffects));

    % histograms
    for i = 1:size(sig2mat,1)
        figure;histogram(sig2mat_perm(i,8,2:end));title(sprintf('Height, %s Model: %s (%.0i permutations)',RandomEffects,RandomEffects{:}(i),nperms));
        figname = sprintf('%s/hist_%s_%s_%s_%.0iperms.pdf',outpath,fname_design(12:end), RandomEffects,RandomEffects{:}(i),nperms);
        saveas(gcf, figname);
    end

    if RandomEffects == "FADSE"
        sig2mat_AD = sig2mat_perm(2,:,:) + sig2mat_perm(3,:,:);
        figure;histogram(sig2mat_AD(:,8,2:end));title('FADSE Model: A+D');
    end

    for i=1:size(fstem_imaging,2)
        load(results_file{i});
        disp(fstem_imaging{i});
        disp(sig2mat);
    
    end

    for i=1:size(fstem_imaging,2)
        load(results_file{i});
        disp(size(iid));
    
    end

end