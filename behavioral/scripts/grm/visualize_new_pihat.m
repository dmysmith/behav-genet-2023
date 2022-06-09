% test script for looking at GRM matrix
filename = '/space/gwas-syn2/1/data/GWAS/ABCD/genotype/release5.0/genotype_QCed/ABCD_20220428.updated.nodups.curated.rel_check.genome';
delimiterIn = ' ';
headerlinesIn = 1;
A = importdata(filename,delimiterIn,headerlinesIn);
colnames = A.textdata(1,:);

% extract vars
IID1 = A.textdata(2:end,2);
IID2 = A.textdata(2:end,4);
PI_HAT = A.data(:,4);

% load iid list from old pihat file
fname_pihat = '/space/amdale/1/tmp/ABCD_cache/abcd-sync/3.0/genomics/ABCD_rel3.0_pihat.mat'; % Should use updated version
load(fname_pihat);
pihatmat(~isfinite(pihatmat)) = 0;


% tmp = [IID1 IID2 num2cell(PI_HAT)]; 
tmp = table(IID1, IID2, PI_HAT); 

% slow method -- loops through all iids
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

% alternate method
pihatmat_new = nan(length(iid_list), length(iid_list));
for i = 1:height(tmp)
    ind1 = strcmp(tmp{i,1},iid_list(:)); 
    ind2 = strcmp(tmp{i,2},iid_list(:));
    pihatmat_new(ind1,ind2) = tmp.PI_HAT(i);
    pihatmat_new(ind2,ind1) = tmp.PI_HAT(i);
end

pihatmat_new(~isfinite(pihatmat_new)) = 0;

% old pihat matrix
addpath(genpath('/home/d9smith/.matlab'));
figure(6); imagesc(pihatmat,0.5*[-1 1]); colormap(blueblackred); axis equal tight; colorbar; xlabel('Subject #'); ylabel('Subject #'); title('Pihat Matrix in ABCD');
export_fig(gcf, '/home/d9smith/tmp/old_pihats.png');

% new pihat matrix
figure(7); imagesc(pihatmat_new,0.5*[-1 1]); colormap(blueblackred); axis equal tight; colorbar; xlabel('Subject #'); ylabel('Subject #'); title('New Pihat Matrix in ABCD');
export_fig(gcf, '/home/d9smith/tmp/new_pihats.png');

% try making figures 
figure(23); hold on;
imagesc(PI_HAT,0.5*[-1 1]); 
colormap(blueblackred); 
axis equal tight;
xlabel('Subject #'); 
ylabel('Subject #'); 
title('Pihat Matrix in ABCD'); hold off;
saveas(gcf,'/home/d9smith/tmp/figure1.png');




% figure(2); hist(colvec(GRM),linspace(-0.5,1.5,201)); xlim([-0.3 1.1]);ylim([0 5000])

