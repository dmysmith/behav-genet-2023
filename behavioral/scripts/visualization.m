% visualization of variances
% diana smith
% may 2022

% inputs
fname_design = 'designMat1_allcovs';
nperms = 1000;
RandomEffects = "FASE";
inpath = "/space/syn50/1/data/ABCD/d9smith/random_effects/behavioral";
outpath = "/home/d9smith/tmp";
phenotypes = {'nihtbx picvocab','nihtbx flanker','nihtbx pattern','nihtbx picture','nihtbx reading','nihtbx cryst','lmt scr perc correct','anthroheightcalc'};


% define set of colors so they match across models
if RandomEffects == "FADSE"
    colors = [0 0.4470 0.7410;0.8500 0.3250 0.0980;0.9290 0.6940 0.1250;0.4940 0.1840 0.5560;0.4660 0.6740 0.1880];
elseif RandomEffects == "FASE"
    colors = [0 0.4470 0.7410;0.8500 0.3250 0.0980;0.4940 0.1840 0.5560;0.4660 0.6740 0.1880];
end

% load results file
results_file = sprintf("%s/%s/nullWB_%dperms/FEMA_wrapper_output_external_%s.mat",inpath,fname_design,nperms,RandomEffects);
load(results_file);

% calculate percentiles
perc = nan(size(sig2mat));
for i = 1:size(sig2mat,1)
    for j = 1:size(sig2mat,2)
        perc(i,j,1) = prctile(sig2mat_perm(i,j,2:end), 2.5);
        perc(i,j,2) = prctile(sig2mat_perm(i,j,2:end), 97.5);
    end
end

% make figure 
model_series = transpose(sig2mat); 
model_error_low = transpose(perc(:,:,1)); 
model_error_high = transpose(perc(:,:,2));
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
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
% Plot the errorbars
errorbar(x',model_series,model_series-model_error_low,model_error_high-model_series,'k','linestyle','none');
hold off
legend(RandomEffects{:}(:));
xticklabels(phenotypes);
title(sprintf('Model including: %s, %s',fname_design(12:end),RandomEffects));

% save figure
figname = sprintf('%s/sig2mat_%s_%s_%.0iperms.png',outpath,fname_design(12:end), RandomEffects,nperms);
saveas(gcf, figname);

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



