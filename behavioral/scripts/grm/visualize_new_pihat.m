% test script for looking at GRM matrix
filename = '/space/gwas-syn2/1/data/GWAS/ABCD/genotype/release5.0/genotype_QCed/ABCD_20220428.updated.nodups.curated.rel_check.genome';
delimiterIn = ' ';
headerlinesIn = 1;
A = importdata(filename,delimiterIn,headerlinesIn);
colnames = A.textdata(1,:);

% extract vars
IID1 = A.textdata(2:end,2);
IID2 = A.textdata(2:end,4);
PI_HAT = A.textdata(2:end,10);

% try making figures 
figure(1); imagesc(GRM,0.5*[-1 1]); colormap(blueblackred); axis equal tight; colorbar; xlabel('Subject #'); ylabel('Subject #'); title('Genetic Relatedness Matrix in ABCD');
figure(2); hist(colvec(GRM),linspace(-0.5,1.5,201)); xlim([-0.3 1.1]);ylim([0 5000])

