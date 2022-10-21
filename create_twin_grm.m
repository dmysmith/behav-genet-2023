%%%%%%%%%%%%%%%%%%%%%
% Create matrix of assigned GRM
% Diana Smith
% Aug 2022

% Add path to cmig_tools (available from:
% https://github.com/cmig-research-group/cmig_tools )
addpath(genpath('/home/d9smith/github/cmig_tools/cmig_tools_utils/matlab'));

% Path where I downloaded some matlab helper functions -- may not be needed for other users
addpath(genpath('/home/d9smith/.matlab'));

% Load full GRM file
grm_file = '/space/abcd-sync/4.0/genomics/ABCD_rel4.0_grm.mat';
load(grm_file);

% load twin ID list
twin_file = "/home/d9smith/projects/random_effects/behavioral/twinfiles/ABCD_twins_all.txt";
twin = readtable(twin_file);
twin = twin(:,["IID1", "IID2", "twin1_genetic_zygosity"]);
twin.measured_grm = repmat(0.5,size(twin,1),1);

% replace with identity matrix
GRM_old = GRM;
GRM = eye(size(GRM));

% create assigned GRM for all participants
GRM_assigned = nan(size(GRM));

for grmi = 1:size(GRM_assigned,1)
    GRM_assigned(grmi,grmi) = 1;
end

% create "notwin" matrix omitting all twin site pairs,
% just to compare distribution of GRM values
GRM_notwin = GRM_old;

% loop through twin pairs
for twini = 1:size(twin,1)
    i = find(strcmp(iid_list,twin{twini,1}));
    j = find(strcmp(iid_list,twin{twini,2}));
    if strcmp(twin{twini,3},'Monozygotic') % assign 1 for mz
        GRM(i,j) = 1;
        GRM(j,i) = 1;
        twin{twini,"measured_grm"} = 1; % if we don't have measured GRM, this will assign 1
        GRM_assigned(i,j) = 1;
        GRM_assigned(j,i) = 1;
    elseif strcmp(twin{twini,3},'Dizygotic') % assign 0.5 for dz
        GRM(i,j) = 0.5;
        GRM(j,i) = 0.5;
        twin{twini,"measured_grm"} = 0.5; % if we don't have measured GRM, this will assign 0.5
    end

    % add measured GRM to twin table
    if isfinite(GRM_old(i,j))
        twin{twini,"measured_grm"}= GRM_old(i,j);
        GRM_notwin(i,j) = NaN;
        GRM_notwin(j,i) = NaN;
    end
end

% save twin GRM files
outfile = '/space/syn50/1/data/ABCD/d9smith/random_effects/behavioral/data/twins_assigned_grm.mat';
save(outfile, 'GRM', 'iid_list');

GRM = GRM_assigned;
assigned_GRM_file = '/space/syn50/1/data/ABCD/d9smith/random_effects/behavioral/data/all_discrete_grm.mat'; 
save(assigned_GRM_file, 'GRM','iid_list')

% for paper -- what is the distribution of GRM in the sample without twin site pts?
for grmi = 1:size(GRM_notwin,1)
    GRM_notwin(grmi,grmi) = NaN;
end

disp(sum(GRM_notwin(:) >= 0.9) / 2);