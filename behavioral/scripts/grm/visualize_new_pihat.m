% test script for looking at GRM matrix

% add paths
addpath(genpath('/home/d9smith/github/cmig_tools/cmig_tools_utils/matlab'));
addpath(genpath('/home/d9smith/.matlab'));

generate_hist=0;

filename = '/space/gwas-syn2/1/data/GWAS/ABCD/genotype/release5.0/genotype_QCed/ABCD_20220428.updated.nodups.curated.rel_check.genome';
delimiterIn = ' ';
headerlinesIn = 1;
A = importdata(filename,delimiterIn,headerlinesIn);

% extract vars
colnames = A.textdata(1,:);
IID1 = A.textdata(2:end,2);
IID2 = A.textdata(2:end,4);
PI_HAT = A.data(:,4);

% load old pihat file
fname_pihat = '/space/amdale/1/tmp/ABCD_cache/abcd-sync/3.0/genomics/ABCD_rel3.0_pihat.mat'; 
load(fname_pihat);
pihatmat(~isfinite(pihatmat)) = 0;


% table of new pihat values 
tmp = table(IID1, IID2, PI_HAT); 

if generate_hist ==1
    % histogram of pihat values
    figure('visible','off');clf; 
    histogram(tmp.PI_HAT, 300);
    % ylim([0 2000]);
    export_fig(gcf, '/home/d9smith/tmp/hist_new_pihats.png');

    % log scale
    figure('visible','off');clf; 
    histogram(tmp.PI_HAT, 1000);
    set(gca,'YScale','log');
    export_fig(gcf, '/home/d9smith/tmp/hist_new_pihats_logscale.png');
end


% populate pihatmat_new with values from tmp.PI_HAT 
pihatmat_new = nan(length(iid_list), length(iid_list));
for i = 1:height(tmp)
    ind1 = strcmp(tmp{i,1},iid_list(:)); 
    ind2 = strcmp(tmp{i,2},iid_list(:));
    pihatmat_new(ind1,ind2) = tmp.PI_HAT(i);
    pihatmat_new(ind2,ind1) = tmp.PI_HAT(i);
end

% everyone should be 1 with themselves
for i = 1:length(iid_list)
    pihatmat_new(i,i) = 1;
end

% save new pihat estimate
outdir = '/space/syn50/1/data/ABCD/d9smith/grm_files';
outfile = matfile(sprintf('%s/%s',outdir,'pihat_20220611.mat'));
outfile.pihatmat = pihatmat_new;

% change all NAs to 0 - is this right?
pihatmat_new(~isfinite(pihatmat_new)) = 0;

if generate_hist == 1
    % plot old pihat matrix
    figure('visible','off');clf; imagesc(pihatmat,0.5*[-1 1]); colormap(blueblackred); axis equal tight; colorbar; xlabel('Subject #'); ylabel('Subject #'); title('Pihat Matrix in ABCD');
    export_fig(gcf, '/home/d9smith/tmp/old_pihats.png');

    % plot new pihat matrix
    figure('visible','off'); clf; imagesc(pihatmat_new,0.5*[-1 1]); colormap(blueblackred); axis equal tight; colorbar; xlabel('Subject #'); ylabel('Subject #'); title('New Pihat Matrix in ABCD');
    export_fig(gcf, '/home/d9smith/tmp/new_pihats.png');

    % what is different between old and new?
    N = 5000;
    idx1 = randperm(length(iid_list),N);
    idx2 = randperm(length(iid_list),N); 

    rand_pihat = [];
    for i = 1:N 
        rand_pihat(i,1) = pihatmat(idx1(i),idx2(i));
        rand_pihat(i,2) = pihatmat_new(idx1(i),idx2(i));
    end

    % plot change from old to new
    figure('visible','off');clf;
    coordLineStyle='k.';
    boxplot(rand_pihat,'Symbol',coordLineStyle); 
    hold on;
    parallelcoords(rand_pihat, 'Color', 0.7*[1 1 1],'LineStyle', '-',...
    'Marker', '.', 'MarkerSize', 10);
    title('change from old to new');
    hold off;
    export_fig(gcf, '/home/d9smith/tmp/pihat_change.png');
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Diana's old code
if 0
    figure(23); hold on;
    imagesc(PI_HAT,0.5*[-1 1]); 
    colormap(blueblackred); 
    axis equal tight;
    xlabel('Subject #'); 
    ylabel('Subject #'); 
    title('Pihat Matrix in ABCD'); hold off;
    saveas(gcf,'/home/d9smith/tmp/figure1.png');
end

% slow method of populating pihatmat_new -- loops through all iids
if 0
    pihatmat_new = nan(length(iid_list), length(iid_list));
    for i=1:length(iid_list)
        for j=1:length(iid_list)
            is_match = strcmp(tmp{:,1},iid_list(i)) & strcmp(tmp{:,2},iid_list(j));
            if(tmp.PI_HAT(is_match))
                pihatmat_new(i,j) = tmp.PI_HAT(is_match);
                pihatmat_new(j,i) = tmp.PI_HAT(is_match); % should go in both directions
            end
        end
    end
end

for i = 1:N
    if rand_pihat(i,1) > .2
        disp(i);
    end
end