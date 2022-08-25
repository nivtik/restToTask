%% Train and save a new model
% All codes are based on Tavor et al., (2016),PMC6309730 and were slightly modified for
% the purposes of the current study
%% Define general variables

% datadir contains rs-fMRI folders for all participants
datadir='/Volumes/IT_RawData/HCP/rfMRI/'; % Path to preprocessed rs-fMRI data folders
outdir='extras/example_output_directory';% Output directory (string)
helper_path='extras/helpers'; % Helpers folder provided in the repository
[cifti,BM] = open_wbfile([helper_path '/example.dtseries.nii']);

% Read test subject lists
test_set=textread(['extras/test_set.txt'],'%s'); %list of test subject numbers
%% Load task contrasts
% One .dtseries file for each task contrast, subjects should be concatenated
% in the same order as in the "subjects" variable

taskdir='Path/to/task/contrasts';
task_list={'Language','WM'};
contrast_list={'MATH_STORY','2BK'};

disp('Loading task contrasts')
test_tasks=zeros(91282,length(test_set),length(contrast_list));

for i=1:length(contrast_list)
    disp(['loading test_set contrast ' task_list{i} '_' contrast_list{i}])
    cifti=open_wbfile([taskdir  '/test_set/AllSubjects_' task_list{i} '_' contrast_list{i} '.dtseries.nii']);
    test_tasks(:,:,i) = cifti.cdata;
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
roi = prod(abs(test_tasks) > 0,2)>0;
ctx = intersect(ctx,find(roi>0));
ctx = intersect(ctx,find(Smask>0));
%% Load group-ICA components

disp('Loading group ICA components')
ICA_cifti=open_wbfile([helper_path '/ica_both_hemis_45_comp.dtseries.nii']);
ICA=ICA_cifti.cdata;

% Determine number of ICA com coponents
numIC=size(ICA,2);

%% Dual regression (creating individaul-level features)
disp('Perfroming DR for test_set');

% HCP session structure (4 resting-state sessions)
root_sessions = {'1' 'LR';'1' 'RL';'2' 'LR';'2' 'RL'};

DR_test_all=zeros(91282,numIC,length(test_set));

for k=1:length(test_set)
    disp(['DR for test subject ' num2str(k)])
    DR_test=DR_both_hemis_single_sub(ICA,datadir,outdir,test_set{k},root_sessions,BM,'save');
    DR_test_all(:,:,k)=DR_test;
end


%% Prediction

for t=1:size(test_tasks,3)
    task_outdir=[outdir '/Prediction_results/' task_list{t} '_' contrast_list{t}];
    
    % load model
    betas=load([outdir '/Prediction_models/' task_list{t} '_' contrast_list{t} '.mat']).betas;
    disp(['Predicting task: ' task_list{t} '_' contrast_list{t}])
    if ~exist(task_outdir)
        mkdir(task_outdir)
    end
    
    % Predict task brain activations
    all_preds=zeros(91282,length(test_set));
    for s=1:length(test_set)
        disp(['Test set subject ' num2str(s)])
        disp('Predict task activity')
        pred=predict_subject_map(DR_test_all(:,:,s),betas,S,task_outdir,helper_path,'dont save');
        all_preds(:,s)=pred;
    end
    
    % Save all predicted maps
    disp(['Saving all preds: ' task_list{t} '_' contrast_list{t}])
    cifti_for_save=cifti;
    cifti_for_save.cdata=all_preds;
    ciftisave(cifti_for_save,[task_outdir '/all_preds.dtseries.nii'])
    
    % Calculate and save predition success measures
    disp('Calculate results')
    [stats,corrmat,null_distribution]=model_stats(test_tasks(:,:,t),all_preds,ctx,1000);
    save([task_outdir '/stats'],'stats')
    save([task_outdir '/CM'],'corrmat')
    save([task_outdir '/null_distribution'])
end