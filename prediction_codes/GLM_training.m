function betas=GLM_training(AllFeatures,all_tasks,outdir,S,helper_path,task_name,save_betas)
% GLM training 

% Created by Saad Jbabdi (Tavor et al., 2016, PMC6309730). Based on Smith
% et al. 2014, PMC4289914.

% Modified to run as a matlab function for usage in: Tik, N, Gal,S, Bernstein-Eliav, M, Tavor,
% I. Towards a generalized AI framework for predicting task-evoked brain 
% activity from resting-state connectivity (2022)

%Input arguments:

% AllFeatures: Individual features (output of dual regression).
% outdir: Output directory
% subjects: Subject list
% helper_path: Path to helper directory
% S: Spatial filters
% task_name: For betas file prefix
% save_betas: Set to "save" in order to save output dtseries.nii file


% Outputs:
% betas: trained model coefficients

addpath(helper_path)

% normalise features
AllFeatures(1:59412,:,:) = normalise(AllFeatures(1:59412,:,:));
AllFeatures(59413:91282,:,:) = normalise(AllFeatures(59413:91282,:,:));
AllFeatures = permute(AllFeatures,[1 3 2]);

task=all_tasks;

% Run training 
unix(['mkdir -p ' outdir '/Betas']);

% start learning
disp('--> start Learning');
featIdx=1:size(AllFeatures,3);
AllFeatures = (normalise(squeeze(AllFeatures(:,:,featIdx))));

% do the GLM
betas=zeros(size(AllFeatures,3)+1,size(S,2));

disp('computing GLM for each parcel')
for j=1:size(S,2)
    disp(j)
    ind = S(:,j)>0;
    y = task(ind,:);
    y=reshape(y, size(y,1)*size(y,2),1);
    
    M=demean(AllFeatures(ind,:,:));
    M=reshape(M,[size(M,1)*size(M,2),size(M,3)]);
    M=[ones(size(M,1),1) M];
    
    betas(:,j) = pinv(M)*y;
end

if strcmp(save_betas,'save')
    save([outdir '/' task_name '.mat'],'betas');
end

end
