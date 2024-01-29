%% Plot figure 6

addpath("../../../shared_functions")

virmen_cell = importdata("../../../shared_files/virmen_cell_5.mat");
tbt_cell = importdata("../../../shared_files/tbt_cell_5.mat");
all_popboot_mean = importdata("all_popboot_means.mat");
all_overall_n_means = importdata("all_overall_n_means.mat");
z_cell_CNN = importdata("../../../shared_files/z_cell_CNN_5_v2.mat");

all_corrs_results = importdata("all_corrs_results_CNN_28072023.mat");

% x_binned_all = importdata("x_binned_all.mat");
% x_binned_all_r = importdata("x_binned_all_r.mat");

accuracies_cell_bvn = importdata("accuracies_cell_neural_50_lr_sep.mat");
accuracies_cell_shuff = importdata("accuracies_cell_neural_50_shuff.mat");

sig_results_all = importdata("sig_results_all.mat");
sig_summary_cell = importdata("sig_summary_cell.mat");
all_full_results_both = importdata("all_full_results_both.mat");
full_all_params_means = importdata("full_all_params_means.mat");
sig_check_all = importdata("sig_check_all.mat");

initialise_params;
linearise_x = true;
bootsamps = 100;
nbins = 50;
xnum = 6;

% [bootstrap_means, centres] = calc_bootstrap_from_raw_subsample_norm(z_cell_CNN{1,1}, virmen_cell{1,1}, tbt_cell{1,1}, linearise_x, nbins,bootsamps, xnum,sub_sample,normalise_z);
centres = importdata("centres_50.mat");
centres = centres*0.74;

zdata = z_cell_CNN{2,2};
virmen_data = virmen_cell{2,2};
tbt_details = tbt_cell{2,2};
sig_results = []; % Why am I not using the sig results??
% sig_results = sig_results_all{2,2};
plot_res = true;
kept_types = [1,3];

% and plot

[wilc_res,h_boots_n] = plot_pop_mean_results_06122022(all_popboot_mean,all_overall_n_means,centres);

plot_correlation_results_06122022(all_corrs_results,centres);

normalise_z = false; % Don't need to normalise as normalised between 0 and 1 anyway
sub_sample = false; % Not subsample, but min trials
[all_sorted_means,all_sorted_means_joined,I_store,cur_neurons_set] = sort_neurons_plot_07122022(zdata,virmen_data,tbt_details,sig_results,linearise_x,nbins,bootsamps,xnum,sub_sample,normalise_z,plot_res,kept_types);

% plot_many_bmi_classifier_results_27022023(x_binned_all,x_binned_all_r,virmen_cell)
[h_boots_svm] = plot_many_svm_bvn_neural_classifier_results(accuracies_cell_bvn,accuracies_cell_shuff,centres);

% Not using examples at the moment
num_examples = 2;
normalise_z = false;
[bootstrap_means, centres] = calc_bootstrap_from_raw_subsample_norm(z_cell_CNN{1,2}, virmen_cell{1,2}, tbt_cell{1,2}, linearise_x, nbins,bootsamps, xnum,sub_sample,normalise_z);
ex_neu = [58,66;164,18;36,110;34,134;17,50;10,24];
ex_dir = [1,1;1,2;1,1;1,1;2,2;1,1];
ex_md = [1,2];
[h_boots] = plot_param_change_figure_18052023(sig_results_all,sig_summary_cell,all_full_results_both,full_all_params_means,sig_check_all,bootstrap_means,centres,ex_neu,ex_dir,ex_md);

plot_across_sessions_mean_change_18042023(all_overall_n_means);
initialise_params;
av_num = 10;
run_many_assess_mean_time_change_18042023(virmen_cell,z_cell_CNN,tbt_cell,model_params,av_num);

full_pairs_means_cell = importdata("full_pairs_means_cell_04072023.mat");

[bootstats_center_kh,bootstats_sem_kh] = plot_many_KH_corrs(full_pairs_means_cell);

all_means_cell = importdata("all_means_cell.mat");
all_stds_cell = importdata("all_stds_cell.mat");

[h_boots2] = plot_decoder_dimension_output(all_means_cell,all_stds_cell);

va_corrs = importdata("va_dec_corrs.mat");
yaw_corrs = importdata("yaw_dec_corrs.mat");
[h_boots_corrs] = plot_weights_changes_corrs(va_corrs,yaw_corrs);