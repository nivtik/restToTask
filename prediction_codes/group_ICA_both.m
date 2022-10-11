function ica_both=group_ICA_both(PCA,cifti,outdir,numIC,save_cifti)

%  Run group ICA on group_pca output (Training set only!)
% 
% Created by Saad Jbabdi (Tavor et al., 2016, PMC6309730). 

% Modified to run as a matlab function for usage in: Tik, N, Bernstein-Eliav, M, Gal,S, Tavor,
% I. Generalizing prediction of task-evoked brain activity across datasets and populations (2022)

% Input arguments:
% PCA: Group PCA results. Could be either path to group PCA results, or
% actual components as a 2D matrix.
% outdir: output directory
% Num_IC: Number of output ICA components (integer)
% save_cifti:  set to "save" in order to save output dtseries.nii file


% Outputs:
% ica_both.dtseries.nii: ICA components calculated on both hemisphres


% set default values for input arguments
    if ~exist('PCA')
        error('Group PCA results not specified');
    end

    if ~exist('outdir')
        outdir='.';
    end
    
    
unix(['mkdir -p ' outdir]);
% Read group PCA results and split Left/Right hemispheres

if isnumeric(PCA)
    Both=double(PCA);
else
cifti=open_wbfile(PCA);
Both   = double(cifti.cdata);
end
%Both=PCA;
ica_both = fastica(Both','approach', 'symm', 'g', 'tanh','lastEig',numIC,'numOfIC', numIC);

% flip sign
ica_both = ica_both .* repmat(sign(sum(sign(ica_both.*(abs(ica_both)>2)),2)),1,size(ica_both,2));

% save
if strcmp(save_cifti,'save')
    dt=cifti;
    dt.cdata=ica_both';
    ciftisave(dt,[outdir '/ica_both_hemis_' num2str(numIC) '_comp.dtseries.nii']);
end
end


