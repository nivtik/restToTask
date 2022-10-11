function DR=DR_both_hemis_single_sub(ICA,datadir,outdir,subj,sessions,BM,save_cifti)

% Run Dual regression to get individual FC components from group ICA

% Created by Saad Jbabdi (Tavor et al., 2016, PMC6309730)

% Modified to run as a matlab function for usage in: Tik, N, Bernstein-Eliav, M, Gal,S, Tavor,
% I. Generalizing prediction of task-evoked brain activity across datasets and populations (2022)


% Input arguments:
% ICA: Group_ICA results, either path to cifti image or data itself
% datadir: directory containing all subjects' resting-state folders
% outdir: output directory
% subj: subject identifier

% Output:
% DR: individual features for each subject

% set default values for input arguments
if ~exist('ICA')
    error('path to group ICA results not defined')
end

if ~exist('datadir')
    error('datadir not defined')
end

if ~exist('outdir')
    outdir=datadir;
end


helper_path='/Volumes/homes/nivtik/Research/Scripts/MATLAB/Predicion_functions/helpers';

unix(['mkdir -p ' outdir]);
addpath(helper_path)

if isempty(subj)
    error('no subjects found')
end


% set defualt value for sessions (HCP data has 4 rs-fMRI sessions)
if ~exist('sessions')
    sessions = {'1' 'LR';'1' 'RL';'2' 'LR';'2' 'RL'};
end

% Dual regression
if ischar(ICA) % if ICA results are saved as a cifti image
    [dt,BM]=open_wbfile(ICA);
    N = size(dt,2);
    G = zeros(91282,N);
    G(BM{1}.DataIndices,1:N) = dt.cdata(BM{1}.DataIndices,:);
    G(BM{2}.DataIndices,1:N) = dt.cdata(BM{2}.DataIndices,:);
    
else % if ICA results are saved as a matrix
    dt=ICA;
    if ~exist('BM')
        [~,BM]=open_wbfile('/Volumes/HCP/HCP_WB_Tutorial_1.0/CP10101_HCP_Pilot-1.fMRI.dtseries.nii');
    end
    N = size(dt,2);
    G = zeros(91282,N);
    G(BM{1}.DataIndices,1:N) = dt(BM{1}.DataIndices,:);
    G(BM{2}.DataIndices,1:N) = dt(BM{2}.DataIndices,:);
end


Hemis = zeros(91282,N);
Hemis(BM{1}.DataIndices,1:N) = 1;
Hemis(BM{2}.DataIndices,1:N) = 1;

pinvG = pinv(G);

disp(subj);

data = [];
disp(' read and demean data');
for sess = 1:size(sessions,1)
    disp(['session ' num2str(sess)]);
    a=sessions{sess,1};b=sessions{sess,2};
    
    % Edit if not using HCP folder structure
    subjdir=[datadir '/' subj '/MNINonLinear/Results' ];
    fname=[subjdir '/rfMRI_REST' a '_' b '/rfMRI_REST' a '_' b '_Atlas_MSMAll_hp2000_clean.dtseries.nii'];
    cifti=open_wbfile(deblank(fname));
    cifti.cdata = variance_normalise(double(cifti.cdata)')';
    data=[data detrend(double(cifti.cdata)')'];
end
% DR - Step 1 (get individual time series)
disp('DR - step 1');
T = pinvG*data;
% DR - Step 2 (get individual spatial maps)
disp('DR - step 2');
[cope,varcope,stats] = fsl_glm(T',data');
DR=stats.t' .* Hemis;

if strcmp(save_cifti,'save')
    % saving subjectwise maps
    DR_dir=[outdir '/DR'];
    if ~exist(DR_dir)
        mkdir(DR_dir)
    end
    oname=[DR_dir '/' subj '_DR2_nosmoothing.dtseries.nii'];
    cifti.cdata = stats.t' .* Hemis;
    ciftisave(cifti,oname);
end
end




