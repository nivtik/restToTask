function pred=predict_subject_map(features,model,S,out_path,helper_path,save_preds)
% inputs - FF file, betas, out_dir, out_name

% Created by Saad Jbabdi (Tavor et al., 2016, PMC6309730). Based on Smith
% et al. 2014, PMC4289914.

% Modified to run as a matlab function for usage in: Tik, N, Gal,S, Bernstein-Eliav, M, Tavor,
% I. Towards a generalized AI framework for predicting task-evoked brain 
% activity from resting-state connectivity (2022)

% Input arguments:
% features: test-set individual features (output of dual regression).
% model: trained model coefficients(betas)
% S: Spatial filters
% out_path: Output path
% helper_path: Path to helper directory
% save_preds: Set to "save" in order to save output dtseries.nii file


% Outputs:
% betas: trained model coefficients

ctx = 1:59412;
subctx = setdiff((1:91282)',ctx);


features(ctx,:) = normalise(features(ctx,:));
features(subctx,:) = normalise(features(subctx,:));
features = normalise(features);

    
if isa(model, 'double')
    pred=zeros(91282,1);
    
    for j=1:size(S,2)
        ind = S(:,j)>0;
        M=demean(features(ind,:));
        M=[ones(size(M,1),1) M];
        pred(ind) = M*model(:,j);
    end
    
elseif isa(model, 'cell')
    pred=zeros(91282,1);
    for j=1:size(S,2)
        ind = S(:,j)>0;
        forest = model{j};
        pred(ind) = predict(forest, features(ind,:));
    end
end

if strcmp(save_preds,'save')
    if ~exist(out_path)
        mkdir(out_path)
    end
    cifti_save = open_wbfile([helper_path '/example.dtseries.nii']);
    cifti_save.cdata = pred;
    ciftisave(cifti_save, [out_path '.dtseries.nii'])
end
end
