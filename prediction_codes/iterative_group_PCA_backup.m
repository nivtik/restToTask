function PCA=iterative_group_PCA(datadir,outdir,subjlist,PCAnew,PCAkeep,sessions, fmri_name, save_cifti)

% Group Incremental PCA (For training set only!)
%  Run PCA on training set data
%
% S.Jbabdi 04/2016

% Modified for usage in: Tik, N, Gal,S, Bernstein-Eliav, M, Tavor,
% I. Towards a generalized AI framework for predicting task-evoked brain 
% activity from resting-state connectivity (2022)

% Input arguments:
% datadir: directory containing all subjects' resting-state folders
% outdir: output directory
% subjlist: a text file containing subject identifiers
% PCAnew: number of components each new iteration contributes
% PCAkeep: number of components to keep
% save_cifti:  set to "save" in order to save output dtseries.nii file

% Output:
% Group_PCA_1000.dtseries.nii
% PCA = group PCA components

% set default values for input arguments
if ~exist('datadir')
    error('datadir not defined')
end

if ~exist('outdir')
    outdir=datadir;
end

if ~exist('subjlist')
    subjects=fdir(datadir,'dir');
end

if ~exist('PCAnew')
    PCAnew=1200;
end

if ~exist('PCAkeep')
    PCAkeep=1000;
end
if ~exist('trim')
    PCAkeep=false;
end

helper_path='/Volumes/homes/nivtik/Research/Scripts/MATLAB/Predicion_functions/helpers';

unix(['mkdir -p ' outdir]);
addpath(helper_path)



% read subject idenifiers from text file (disable if subjects received as

if iscell(subjlist)
    subjects = subjlist;
else
    subjects = textread(subjlist, '%s');
end

if isempty(subjects)
    error('no subjects found')
end

% Loop over sessions and subjects
W= [] ;

for s=1:length(subjects)
    subj=subjects{s};
    disp(subj);
    
    % If not using HCP folder structure, chagne to match your own
    subjdir=[datadir '/' subj '/MNINonLinear/Results/'];
    fname=[subjdir fmri_name  '/' fmri_name '_Atlas_MSMAll_hp2000_clean.dtseries.nii'];
    
    % read and demean data
    disp('read data');
    [cifti,BM]=open_wbfile(deblank(fname));
    grot=demean(double(cifti.cdata)'); clear cifti.cdata;
    
    % noise variance normalisation
    grot = variance_normalise(grot);
 
    % concat
    W=[W; demean(grot)]; clear grot;
    % PCA reduce W to PCAnew eigenvectors
    disp(['do PCA ' num2str(size(W,1)) 'x' num2str(size(W,2))]);
    [uu,dd]=eigs(W*W',min(PCAnew,size(W,1)-1));  W=uu'*W; clear uu;
end


data=W(1:min(PCAkeep,size(W,1)),:)';

% Output PCA components
dt=open_wbfile([helper_path '/example.dtseries.nii']);
dt.cdata=data;
PCA=data;

% Save PCA components as cifti
if strcmp(save_cifti,'save')
    % Save group PCA results
    ciftisave(dt,[outdir '/Group_PCA_' num2str(PCAkeep) '.dtseries.nii']);
end

end





