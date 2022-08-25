function PCA=group_PCA(datadir,outdir,subjlist,PCAnew,PCAkeep, sessions, fmri_name, helper_path, save_cifti)
% Group Incremental PCA (For training set only!)
%  Run PCA on training set data

% Created by Saad Jbabdi (Tavor et al., 2016, PMC6309730). Based on Smith
% et al. 2014, PMC4289914.

% Modified to run as a matlab function for usage in: Tik, N, Gal,S, Bernstein-Eliav, M, Tavor,
% I. Towards a generalized AI framework for predicting task-evoked brain 
% activity from resting-state connectivity (2022)

% Input arguments:
% datadir: directory containing all subjects' resting-state folders
% outdir: output directory
% subjlist: a text file containing subject identifiers
% PCAnew: number of components each new iteration contributes
% PCAkeep: number of components to keep
% sessions: rs-fMRI session prefixes
% save_cifti:  set to "save" in order to save output dtseries.nii file

% Output:
% Group_PCA_1000.dtseries.nii
% PCA = group PCA components

% set default values for input arguments
if ~exist('datadir')
    error('datadir not defined')
end
if ~exist('save_cifti')
    save_cifti='save';
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


unix(['mkdir -p ' outdir]);


% read subject idenifiers from text file
if exist('subjlist')
    if ischar(subjlist)
        subjects = textread(subjlist, '%s');
    else
        subjects = subjlist;
end

if isempty(subjects)
    error('no subjects found')
end

% Loop over sessions and subjects
W= [] ;
for sess = 1:size(sessions,1)
    a=sessions{sess,1};b=sessions{sess,2};
    
    for s=1:length(subjects)
        subj=subjects{s};
        disp(subj);
        
        % Change below if not working with HCP folder structure
        subjdir=[datadir '/' subj '/MNINonLinear/Results/' ];
        fname=[subjdir '/' fmri_name a '_' b '/' fmri_name a '_' b '_Atlas_MSMAll_hp2000_clean.dtseries.nii'];
        
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
    
end
data=W(1:PCAkeep,:)';

% Output PCA components
PCA=data;

if strcmp(save_cifti,'save')
    % Save group PCA results
    dt=open_wbfile([helper_path '/example.dtseries.nii']);
    dt.cdata=data;
    PCA_cifti=dt;
    ciftisave(dt,[outdir '/Group_PCA_' num2str(PCAkeep) '.dtseries.nii']);
end

end