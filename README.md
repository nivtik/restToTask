The code in this repository is associated with the following paper: Tik, N, Bernstein-Eliav, M, Gal,S, Tavor,I. Generalizing prediction of task-evoked brain activity across datasets and populations (2022)

All feature extraction, model training and prediction codes are based on the codes released with: Tavor, I., Jones, O. P., Mars, R. B., Smith, S. M., Behrens, T. E., & Jbabdi, S. (2016). Task-free MRI predicts individual differences in brain activity during task performance. Science, 352(6282), 216-220. Original codes could be found in https://git.fmrib.ox.ac.uk/saad/ActPred

The following software are required:
  1. MATLAB (https://www.mathworks.com/products/matlab.html)
  2. FSL (https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/)
  3. Connectome Workbench (https://www.humanconnectome.org/software/get-connectome-workbench)
  4. FastICA toolbox (http://research.ics.aalto.fi/ica/fastica)

All the main scripts and additional code are provied here.
The main dataset is provided by the Human Connectome Project via https://db.humanconnectome.org.

Provided resources:
   - Seed group-ICA map that can be used for trainng new models: extras/helpers/ica_both_hemis_45_comp.dtseries.nii
   - Trained models (beta coefficients) for 18 HCP task contrasts (740 training participants,1200 timepoints each): extras/HCP_YA_trained_models
     
  1. Calculate group-ICA components on new data:
  
    - Main script: genertate_group_components.m
    - Additional code: group_PCA.m, group_ICA_both.m

   2. Train a new model:
   
    - Main script: model_training.m
    - Additional code: DR_both_hemis_single_sub.m, GLM_training.m

   3. Prediction:
   
    - Main script: prediction_with_pretrained_model.m
    - Additional code: DR_both_hemis_single_sub.m, predict_subject_map.m, model_stats.m

