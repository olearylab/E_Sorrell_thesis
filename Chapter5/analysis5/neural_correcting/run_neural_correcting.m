%% Run neural correlations with error correcting movements

% load in virmen and tbt_cell and z_cell_CNN
% get z_binned_cell using e.g.
% Used normalise_z = true; nbins = 50;
% [z_binned_cell,error_cell,z_error_corrs_cell,z_error_peak_corrs_cell] = run_many_compare_neural_to_error(z_cell,virmen_cell,tbt_cell,model_params,normalise_z,nbins);

% Also get error_cell and x_binned_cell using e.g.
% [error_cell,x_cell,ro_all,sign_matching,h_boots] = run_many_check_BMI_correct_errors(virmen_cell,tbt_cell,nbins,false)

[z_correction_corrs_cell,z_correction_corrs_p_cell] = run_many_neuron_correction_corrs(z_binned_cell,x_binned_cell,error_cell,tbt_cell);

% To run cross validated decoders. load in train_mean_cell and Rfixed_cell
initialise_params;
[all_res_cell,ball_weights_cell,bmi_weights_cell] = run_many_decoding_check_pixels(virmen_cell,train_mean_cell,Rfixed_cell,tbt_cell,model_params);

