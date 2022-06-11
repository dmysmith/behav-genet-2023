% note from Diana -- need to add anders function to path
addpath(genpath('/home/dale/matlab/MOSTest'));
addpath(genpath('/home/dale/matlab/utils'));
addpath(genpath('/home/d9smith/.matlab'));
addpath(genpath('/home/d9smith/github/cmig_tools/cmig_tools_utils'));

try
   canUseGPU = parallel.gpu.GPUDevice.isAvailable;
catch ME
   canUseGPU = false;
end

bfile = '/space/gwas-syn2/1/data/GWAS/ABCD/RUCDR_batch/MERGE_202204/ABCD_20220428.updated.nodups.curated';
fname_bim = sprintf('%s.bim',bfile);
%fname_bim = '/space/gwas-syn2/1/data/GWAS/ABCD/RUCDR_batch/MERGE_202204/ABCD_20220428.updated.nodups.curated.bim'; % Chun's updated genotypes
fname_fam = strrep(fname_bim,'.bim','.fam');
%fname_fam = '/space/syn09/1/data/rloughna/ABCD/MOSTest/ABCD_genotypes/HRC/merged_ABCD.fam';
fileID = fopen(fname_bim);
bim_file = textscan(fileID,'%s %s %s %s %s %s');
fclose(fileID);
snps=length(bim_file{1});
fileID = fopen(fname_fam);
fam_file = textscan(fileID,'%s %s %s %s %s %s');
fclose(fileID);
nsubj=length(fam_file{1});
PIDs_geno = fam_file{1};
FIDs_geno = fam_file{2};
fprintf('%s: %i snps and %i subjects detected in bimfile %s and famfile %s\n', bfile, snps, nsubj, fname_bim, fname_fam);

% geno_int8 = PlinkRead_binary2(nsubj, 1:snps, bfile);
geno_int8 = PlinkRead_binary2_amd(nsubj, 1:snps, bfile);
geno = nan(size(geno_int8), 'single'); for code = int8([0,1,2]), geno(geno_int8==code) = single(code); end;

mafvec = nanmean(geno,1)/2; mafvec = min(mafvec,1-mafvec);
cratevec = mean(isfinite(geno),1); % Why are these all either 0 or 1? Check if imputation includes low-confidence calls

ivec_snps = find(cratevec>0.95 & mafvec>0.05);
geno_tmp = (geno-nanmean(geno))./nanstd(geno);
geno_tmp(~isfinite(geno_tmp)) = 0;
if canUseGPU
  geno_tmp = gpuArray(geno_tmp);
end

% Compute covariance matrix from z-transformed SNP dose matrix
tic
grm_sum = geno_tmp(:,ivec_snps)*geno_tmp(:,ivec_snps)'; nsum = length(ivec_snps);
toc

GRM = grm_sum/nsum;
[U S] = svd(GRM); s = diag(S);
%n = 10; GRM_corr = cov2corr(GRM - U(:,1:n)*diag(s(1:n))*U(:,1:n)');
n = 10; GRM_corr = cov2corr(GRM - U(:,1:n)*diag(s(1:n))*U(:,1:n)');
GRM = cov2corr(GRM);
figure('visible','off'); clf; imagesc(GRM,0.5*[-1 1]); colormap(blueblackred); axis equal tight; colorbar; xlabel('Subject #'); ylabel('Subject #'); title('Genetic Relatedness Matrix in ABCD');export_fig(gcf, '/home/d9smith/tmp/grm.png');
figure('visible','off'); clf; hist(colvec(GRM),linspace(-0.5,1.5,201)); xlim([-0.3 1.1]);ylim([0 5000]);export_fig(gcf, '/home/d9smith/tmp/fig2.png');
figure('visible','off'); clf; imagesc(GRM_corr,0.5*[-1 1]); colormap(blueblackred); axis equal tight; colorbar; xlabel('Subject #'); ylabel('Subject #'); title('PCA-Corrected Genetic Relatedness Matrix in ABCD');export_fig(gcf, '/home/d9smith/tmp/grm_pca.png');
figure('visible','off'); clf; hist(colvec(GRM_corr),linspace(-0.5,1.5,201)); xlim([-0.3 1.1]);ylim([0 5000]);export_fig(gcf, '/home/d9smith/tmp/fig12.png');

% save out GRM and GRM_corr
outdir = '/space/syn50/1/data/ABCD/d9smith/grm_files';
outfile = matfile(sprintf('%s/%s',outdir,'grm.mat'));
outfile.GRM = GRM;
outfile.GRM_corr = GRM_corr;


% Use svd of z-transformed SNP dose matrix
tic
[U S V] = svd(geno_tmp,'econ'); s = diag(S); % This is kind of slow -- should try on GPUs?
toc
% Look at the covariance of the first few principal components
n = 3; tmp = U(:,1:n)*diag(s(1:n).^2)*U(:,1:n)';
figure('visible','off'); clf; imagesc(tmp,max(abs(tmp(:)))*[-1 1]); colormap(blueblackred); axis equal tight; colorbar; xlabel('Subject #'); ylabel('Subject #');export_fig(gcf, '/home/d9smith/tmp/fig666.png');


% Look at batch effects
tbl_batch = readtable('/space/gwas-syn2/1/data/GWAS/ABCD/genotype/release3.0/ABCD_release3.0_.batch_info.txt');
[dummy IA IB] = intersect(tbl_batch.abcd_id_redcap,PIDs_geno,'stable');
[Axiom_list dummy2 IC_Axiom] = unique(tbl_batch.Axiom_Plate);
[BATCH_list dummy2 IC_BATCH] = unique(tbl_batch.BATCH);

batchvec = NaN(1,nsubj); batchvec(IB) = IC_BATCH(IA);
for bi = 1:length(BATCH_list)
  mv = mean(GRM(:,batchvec==bi));
  fprintf(1,'batch %d (%s): mean(GRM) = %f\n',bi,BATCH_list(bi),mv);
end


% Look at PCs in tabulated data
fname_covars = '/space/amdale/1/tmp/ABCD_cache/abcd-sync/3.0/support_files/design_matrices/ABCD_3.0_design_example1.txt'; % Should be passed in to function, update with covars only, 4.0 data, age & age^2, imputed
tbl_covars = readtable(fname_covars);
[dummy IA IB] = intersect(tbl_covars.src_subject_id,PIDs_geno,'stable');

tmp = table2array(tbl_covars(IA,[7:16]));
GRM_PCs = cov2corr(nancov(tmp'));

figure(21); imagesc(GRM(IB,IB),0.5*[-1 1]); colormap(blueblackred); axis equal tight; colorbar; xlabel('Subject #'); ylabel('Subject #'); title('Genetic Relatedness Matrix in ABCD');
figure(31); imagesc(GRM_PCs,1.0*[-1 1]); colormap(blueblackred); axis equal tight; colorbar; xlabel('Subject #'); ylabel('Subject #'); title('Genetic Relatedness Matrix in ABCD');


% ToDo
%   Look at genotyping batch effects: /space/gwas-syn2/1/data/GWAS/ABCD/genotype/release3.0/ABCD_release3.0_.batch_info.txt

%   Identify subjects by predominant ancestry, look at within vs. across ancestry grm
%   Look at cross-site vs. within site (latter is more likely to have close relatives)

%   Check into why standard grm formula doesn't quite work (results in GRM values much beyond 1.0)



% Examine Chun's pihat estimates 

fname_pihat = '/space/amdale/1/tmp/ABCD_cache/abcd-sync/3.0/genomics/ABCD_rel3.0_pihat.mat'; % Should use updated version
load(fname_pihat);
pihatmat(~isfinite(pihatmat)) = 0;

figure(101); imagesc(pihatmat,0.5*[-1 1]); colormap(blueblackred); axis equal tight; colorbar; xlabel('Subject #'); ylabel('Subject #'); title('Pihat Matrix in ABCD');

[U S] = svd(pihatmat-eye(size(pihatmat)));
s = diag(S); 
n = 3; pihatmat_corr = cov2corr(pihatmat - U(:,1:n)*diag(s(1:n))*U(:,1:n)');
sfigure(111); imagesc(pihatmat_corr,0.5*[-1 1]); colormap(blueblackred); axis equal tight; colorbar; xlabel('Subject #'); ylabel('Subject #'); title('PCA-Corrected Pihat Matrix in ABCD');



% Extras / old junk

dirname_geno_cache = '/space/amdale/1/tmp/GWAS_cache/ABCD';
flist = dir(sprintf('%s/genomat_chunk*',dirname_geno_cache)); flist = {flist.name}; nfiles = length(flist);

try
   canUseGPU = parallel.gpu.GPUDevice.isAvailable;
catch ME
   canUseGPU = false;
end

grm_sum = 0; nsum = 0;
for filei = 1:nfiles
  if mod(filei,10)==0
    fprintf(1,'filei=%d/%d (%s)\n',filei,nfiles,datestr(now));
    GRM = cov2corr(grm_sum/nsum);
%    GRM = grm_sum/nsum;
    sfigure(1); imagesc(GRM,0.5*[-1 1]); colormap(blueblackred); axis equal tight; colorbar; xlabel('Subject #'); ylabel('Subject #'); title('Genetic Relatedness Matrix in ABCD');
    sfigure(2); hist(colvec(GRM),linspace(-0.5,1.5,201));
%    sfigure(3); hist(diag(GRM),linspace(-0.5,1.5,201)); % Should be all ones (subject has to be perfectly correllated with himself)
    [U S] = svd(GRM-eye(size(GRM)));
    s = diag(S);
    n = 3; GRM_corr = cov2corr(GRM - U(:,1:n)*diag(s(1:n))*U(:,1:n)');
    sfigure(11); imagesc(GRM_corr,0.5*[-1 1]); colormap(blueblackred); axis equal tight; colorbar; xlabel('Subject #'); ylabel('Subject #'); title('PCA-Corrected Genetic Relatedness Matrix in ABCD');
    sfigure(12); hist(colvec(GRM_corr),linspace(-0.5,1.5,201));
    drawnow;
  end
  load(sprintf('%s/%s',dirname_geno_cache,flist{filei}));
%  geno = double(geno);
  mafvec = nanmean(geno,1)/2; mafvec = min(mafvec,1-mafvec);
  callvec = mean(isfinite(geno),1); % Why are these all either 0 or 1? Check if imputation includes low-confidence calls
  ivec_snps = find(callvec>0.95 & mafvec>0.05);
%  geno_tmp = (geno-2*mafvec)./sqrt(mafvec.*(1-mafvec));
  geno_tmp = (geno-mean(geno))./std(geno);
  if canUseGPU
    geno_tmp = gpuArray(geno_tmp);
  end
%  geno_tmp(~isfinite(geno_tmp)) = 0;
  grm_sum = grm_sum + geno_tmp(:,ivec_snps)*geno_tmp(:,ivec_snps)'; nsum = nsum + length(ivec_snps);
end

