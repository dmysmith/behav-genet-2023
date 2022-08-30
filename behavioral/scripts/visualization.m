% visualization of variances
% diana smith
% may 2022

% inputs
results_dir = '/space/syn50/1/data/ABCD/d9smith/random_effects/behavioral/results/results_20220829';
modelname = {'model1'; 'model2'; 'model3'; 'model4'};

outpath = "/home/d9smith/projects/random_effects/behavioral/results/plots"; % outpath = "/home/d9smith/tmp";
results_file = strcat(results_dir,'/FEMA_wrapper_output_external_',modelname,'.mat');

titles = {
'FAE Model, twins only at baseline, GRM assumed';
'FAE Model, twins only at baseline, with GRM';
'FAE Model, full sample at baseline, with GRM';
'FATE Model, full sample at baseline, with GRM';
'FASTE Model, full sample at baseline and year 2, with GRM';
'FAE Model, full sample at baseline, discrete zygosity';
'FATE Model, full sample at baseline, discrete zygosity';
'FASTE Model, full sample at baseline and year 2, discrete zygosity';
'FAE Model, twins only at baseline, GRM assumed, residualized for covariates';
'FAE Model, twins only at baseline, with GRM, residualized for covariates';
'FAE Model, full sample at baseline, with GRM, residualized for covariates';
'FATE Model, full sample at baseline, with GRM, residualized for covariates';
'FASTE Model, full sample at baseline and year 2, with GRM, residualized for covariates';
'FASE model, full sample minus twins at baseline and year 2, with GRM';
}


for i=1:size(modelname,1)
    load(results_file{i});
    phenotypes = colnames_imaging;

    % define set of colors so they match across models
    if isequal(RandomEffects{i}, {'F';'A';'D';'S';'E'})
        colors = [0 0.4470 0.7410;0.8500 0.3250 0.0980;0.9290 0.6940 0.1250;0.4940 0.1840 0.5560;0.4660 0.6740 0.1880];
    elseif isequal(RandomEffects{i}, {'F';'A';'S';'E'})
        colors = [0 0.4470 0.7410;0.8500 0.3250 0.0980;0.4940 0.1840 0.5560;0.4660 0.6740 0.1880];
    elseif isequal(RandomEffects{i}, {'F';'A';'E'})
        colors = [0 0.4470 0.7410;0.8500 0.3250 0.0980;0.4660 0.6740 0.1880]; 
    end

    % make figure 
    model_series = transpose(sig2mat(:,:,1)); 
    model_error_low = transpose(sig2mat(:,:,2)); 
    model_error_high = transpose(sig2mat(:,:,3));
    b = bar(model_series,'FaceColor','flat');
    for k = 1:size(sig2mat,1)
        b(k).CData = colors(k,:);
    end
    xlabel('Phenotype');
    ylabel('Percent of Variance');
    hold on
    % Calculate the number of groups and number of bars in each group
    [ngroups,nbars] = size(model_series);
    % Get the x coordinate of the bars
    x = nan(nbars, ngroups);
    for bari = 1:nbars
        x(bari,:) = b(bari).XEndPoints;
    end
    % Plot the errorbars
    errorbar(x',model_series,model_series-model_error_low,model_error_high-model_series,'k','linestyle','none');
    hold off
    legend(RandomEffects{i}(:), 'Location','eastoutside');
    set(gca,'TickLabelInterpreter','none')
    xticklabels(phenotypes);
    xtickangle(45);
    title(sprintf('%s: %s',fstem_imaging{i}, titles{i}));

    % save figure
    figname = sprintf('%s/%s.png',outpath,fstem_imaging{i});
    saveas(gcf, figname);

end

%% old code
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

end