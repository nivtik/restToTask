%% Train and save a new model

% All codes are based on Tavor et al., (2016),PMC6309730 and were slightly modified for
% the purposes of: Tik, N, Bernstein-Eliav, M, Gal,S, Tavor, I. Generalizing 
% prediction of task-evoked brain activity across datasets and populations (2022)

%% Define general variables

% datadir contains rs-fMRI folders for all participants
datadir='/Volumes/IT_RawData/HCP/rfMRI/'; % Path to preprocessed rs-fMRI data folders
outdir='extras/example_output_directory';% Output directory (string)
helper_path='extras/helpers'; % Helpers folder provided in the repository
addpath(helper_path);
[cifti,BM] = open_wbfile([helper_path '/example.dtseries.nii']);

% Read training subject lists
subjects=textread('extras/training_set.txt','%s'); %list of subject numbers
%% Load task contrasts
% One .dtseries file for each task contrast, subjects should be concatenated
% in the same order as in the "subjects" variable

taskdir='Path/to/task/contrasts';
task_list={'Language','WM'};
contrast_list={'MATH_STORY','2BK'};

disp('Loading task contrasts')
training_tasks=zeros(91282,length(subjects),length(contrast_list));

for i=1:length(contrast_list)
    disp(['loading training_set contrast ' task_list{i} '_' contrast_list{i}])
    cifti=open_wbfile([taskdir  '/training_set/AllSubjects_' task_list{i} '_' contrast_list{i} '.dtseries.nii']);
    training_tasks(:,:,i) = cifti.cdata;
end

%% Load spatial Smasks
disp('Load Smasks');
Smask = open_wbfile([helper_path '/ica_both_lowdim_subjectsE.dtseries.nii']);
[m,wta]=max(Smask.cdata,[],2);
wta = wta .* (m>2.1);
S = zeros(size(Smask.cdata));
for i=1:size(Smask.cdata,2)
    S(:,i) = double(wta==i);
end

Smask=sum(S,2);
%% ctx

ctx = [BM{1}.DataIndices(:);BM{2}.DataIndices(:)];
subctx = setdiff((1:91282)',ctx);
Lctx=BM{1}.DataIndices(:);
Rctx=BM{2}.DataIndices(:);

% remove empty grayordinates
disp('removing empty grayordiantes')
roi = prod(abs(training_tasks) > 0,2)>0;
ctx = intersect(ctx,find(roi>0));
ctx = intersect(ctx,find(Smask>0));
%% Load group-ICA components

disp('Loading group ICA components')
ICA_cifti=open_wbfile([helper_path '/ica_both_hemis_45_comp.dtseries.nii']);
ICA=ICA_cifti.cdata;

% Determine number of ICA com coponents
numIC=size(ICA,2);

%% Dual regression (creating individaul-level features)
disp('Perfroming DR for training_set');

% HCP session structure (4 resting-state sessions)
root_sessions = {'1' 'LR';'1' 'RL';'2' 'LR';'2' 'RL'};

DR_training_all=zeros(91282,numIC,length(subjects));

for k=1:length(subjects)
    disp(['DR for training subject ' num2str(k)])
    DR_training=DR_both_hemis_single_sub(ICA,datadir,outdir,subjects{k},root_sessions,BM,'save');
    DR_training_all(:,:,k)=DR_training;
end

%% Model training

for t=1:size(training_tasks,3)
    task_outdir=[outdir '/Prediction_models'];
    
    disp(['Predicting task: ' task_list{t} '_' contrast_list{t}])
    if ~exist(task_outdir)
        mkdir(task_outdir)
    end
    
    % Train the prediction model
    GLM_training(DR_training_all,training_tasks(:,:,t),task_outdir,S,helper_path,[task_list{t} '_' contrast_list{t}],'save');
    
end
