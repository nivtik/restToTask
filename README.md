The code in this repository is associated with the following publication: Tik, N, Gal,S, Bernstein-Eliav, M, Tavor,
% I. Towards a generalized AI framework for predicting task-evoked brain 
% activity from resting-state connectivity (2022)

The following software are required:
  1. MATLAB (https://www.mathworks.com/products/matlab.html)
  2. FSL (https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/)
  3. Connectome Workbench (https://www.humanconnectome.org/software/get-connectome-workbench)
  4. FastICA toolbox (http://research.ics.aalto.fi/ica/fastica)

All the main scripts and additional code are provied here.
The main dataset is provided by the Human Connectome Project via https://db.humanconnectome.org.

Provided resources:
   - Seed group-ICA map that can be used for trainng new models: extras/
   - Trained models (beta coefficients) for 18 HCP task contrasts detailed in the paper
   
Provided code:
  1. Calculate group-ICA components on new data:
    - Main script: genertate_group_components.m
    - Additional code: group_PCA.m, group_ICA_both.m

   2. Train a new model:
    - Main script: model_training.m
    - Additional code: DR_both_hemis_single_sub.m, GLM_training.m

   3. Prediction:
    - Main script: prediction_with_pretrained_model.m
    - Additional code: DR_both_hemis_single_sub.m, predict_subject_map.m, model_stats.m
