%% Run example DFF filt

% Load Data for DW48
% Need virmen, neural and tbt_details

%

% Only 60 training trials
training_trials = 60;

[ztrain, ztest, xtrain, xtest] = split_session_data(zfull, xfull, tbt_details, [1,2], training_trials, true, false);

% raw
initialise_params;
model_params.reg_images = false;
model_params.spatial = false;
model_params.dff = false;
[results_struct] = DW_train_test_offline_func_05_2021(ztrain, ztest, xtrain, xtest, model_params, tbt_details, [1,2], 7, 4);

% dff filtered
model_params.dff = true;
[results_struct_raw] = DW_train_test_offline_func_05_2021(ztrain, ztest, xtrain, xtest, model_params, tbt_details, [1,2], 7, 4);