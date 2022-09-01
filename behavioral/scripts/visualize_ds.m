function b = visualize_ds(results_file, outpath, fstem_imaging, titles, i)
    % visualize_ds visualize FEMA random effects estimates with confidence intervals
    % Saves a file with a bar chart of FEMA random effects estimates.

    load(results_file);
    phenotypes = colnames_imaging;

    % define set of colors so they match across models
    if isequal(RandomEffects{i}, {'F';'A';'D';'S';'E'})
        colors = [0 0.4470 0.7410;0.8500 0.3250 0.0980;0.9290 0.6940 0.1250;0.4940 0.1840 0.5560;0.4660 0.6740 0.1880];
    elseif isequal(RandomEffects{i}, {'F';'A';'S';'E'})
        colors = [0 0.4470 0.7410;0.8500 0.3250 0.0980;0.4940 0.1840 0.5560;0.4660 0.6740 0.1880];
    elseif isequal(RandomEffects{i}, {'F';'A';'E'})
        colors = [0 0.4470 0.7410;0.8500 0.3250 0.0980;0.4660 0.6740 0.1880];
    elseif isequal(RandomEffects{i}, {'F';'A';'T';'E'})
        colors = [0 0.4470 0.7410;0.8500 0.3250 0.0980;0.9290 0.6940 0.1250;0.4660 0.6740 0.1880]; 
    elseif isequal(RandomEffects{i}, {'F';'A';'S';'T';'E'})
        colors = [0 0.4470 0.7410;0.8500 0.3250 0.0980;0.4940 0.1840 0.5560;0.9290 0.6940 0.1250;0.4660 0.6740 0.1880]; 
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
    ylim([0 1]);
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
    title(sprintf('%s: %s',fstem_imaging{i}, titles{i}),'interpreter', 'none');

    % save figure
    figname = sprintf('%s/%s.png',outpath,fstem_imaging{i});
    saveas(gcf, figname);
    fprintf('Saved figure to %s\n', figname);

return