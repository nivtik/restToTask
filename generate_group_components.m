%% Generate group-ICA compontents for prediction
% All codes are based on Tavor et al., (2016),PMC6309730 and were slightly modified for
% the purposes of the current study

% Group-ICA components could be calculated using the training set or an outgroup
% Do not include any test set data!

%% Define general variables

% datadir contains rs-fMRI folders for all participants. 
% group-ICA maps in the extras folder were computed using 100 unrelated
% participants that could be downloaded from Connectome DB
datadir='/Volumes/IT_RawData/HCP/rfMRI/'; 
outdir='./extras/example_output_directory';% Output directory (string)
helper_path='./extras/helpers'; % Helpers folder provided in the repository
addpath(helper_path);
[cifti,BM] = open_wbfile([helper_path '/example.dtseries.nii']);

% Read training subject lists
subjects=textread(['./extras/100_unrelated.txt'],'%s'); %list of test subject numbers
%% Iterative group PCA (for dimensionality reduction)
fmri_name='rfMRI_REST';
sessions={'1' 'LR';'1' 'RL';'2' 'LR';'2' 'RL'}; % For HCP data
PCA=group_PCA(datadir,outdir,subjects,1200,1000,sessions, fmri_name,helper_path, 'save');
%% Group-ICA for generating components

numIC=45; %Number of requested ICA components 
ica_both=group_ICA_both(PCA,cifti,outdir,numIC,'save');