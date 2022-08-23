%%%%%%%%%%%%%%%%%%%%%
% Create dummy GRM for twin analyses
% Diana Smith
% Aug 2022

% Load full GRM file
grm_file = '/space/amdale/1/tmp/ABCD_cache/abcd-sync/4.0/genomics/ABCD_rel4.0_grm.mat';
load(grm_file);

% load twin ID list
twin_file = '/space/syn50/1/data/ABCD/d9smith/random_effects/behavioral/data/pheno/twins_ids_for_grm.txt';
twin = readtable(twin_file);

% replace with identity matrix
GRM_old = GRM;
GRM = eye(size(GRM));

% loop through twin pairs
for twini = 1:size(twin,1)
    i = find(strcmp(iid_list,twin{twini,1}));
    j = find(strcmp(iid_list,twin{twini,2}));
    if strcmp(twin{twini,3},'Monozygotic') % assign 1 for mz
        GRM(i,j) = 1;
        GRM(j,i) = 1;
    elseif strcmp(twin{twini,3},'Dizygotic') % assign 0.5 for dz
        GRM(i,j) = 0.5;
        GRM(j,i) = 0.5;
    end
end

% save new GRM file
outfile = '/space/syn50/1/data/ABCD/d9smith/random_effects/behavioral/data/twins_assigned_grm.mat';
save(outfile, 'GRM', 'iid_list');