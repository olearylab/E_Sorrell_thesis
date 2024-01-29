%% Run neural change full

% Mean population changes
% load in z_cell_CNN, virmen_cell, tbt_cell
initialise_params;
[all_popboot_mean,all_neuron_means,all_overall_n_means] = run_many_mean_population_change(z_cell_CNN,virmen_cell,tbt_cell,model_params,true,50,100,6,true,[1,4;7,10]);

% Run neural svm classifiers
[accuracies_cell] = run_many_svm_bvn_neural_classifier(z_cell,virmen_cell,tbt_cell,50,[1,2,3,4],true,false,true);
% and shuffles
n_shuffles = 100;
[accuracies_cell_shuffle] = run_many_svm_bvn_neural_classifier_shuffle(z_cell,virmen_cell,tbt_cell,50,[1,2,3,4],true,false,true,n_shuffles);

% Run significant tuning check
shuff_limit = 300;
num_shuffles = 1000;
[sig_results_all,sig_check_all] = run_many_sig_tuning_check(z_cell, virmen_cell, tbt_cell, num_shuffles, true, 50, 100, shuff_limit, 6);
[sig_summary_cell] = calc_peak_gain_loss(sig_results_all);

% Run tuning curve changes
sub_sample = false;
[all_proportions,all_proportions_overlap,all_proportions_both,all_full_results,all_full_results_overlap,all_full_results_both,all_mean_curve_params,full_all_params] = assess_tc_param_changes_many_21022023(z_cell,virmen_cell,tbt_cell,true,50,100,[2.5,97.5],6,sig_check_all,sub_sample);

% Run noise correlations
sub_sample = false;
[full_pairs_means_cell] = run_many_KH_correlations_04072023(z_cell,virmen_cell,tbt_cell, model_params, 50, true, sub_sample, true);

% check projection through decoder dimenions
% Need n_W_online and offline_cell. These are the weights within neurons
% calculated from the pixel weights. I'm pretty sure I used this, but I
% also have function for the projections through the entire weights.
[all_means_cell,all_means_cell_ones,all_stds_cell] = run_many_decoder_dimension_check_n(n_W_online_cell,n_W_offline_cell,z_cell,virmen_cell,tbt_cell,model_params);

% Calculate correlations of weights and activity changes
% need stat cell, online weights, offline bmi trial weights, and neuron
% changes in mean activity, and neurons stds (although I think I don't
% actually use their scaling.
[va_corrs,yaw_corrs] = compare_weights_to_changes(Wout_online_cell,bmi_weights_cell,stat_cell,all_overall_n_means,z_stds);