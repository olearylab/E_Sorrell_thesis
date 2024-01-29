%% Plot neural correcting

addpath("../../../shared_functions")

virmen_cell = importdata("../../../shared_files/virmen_cell_5.mat");
tbt_cell = importdata("../../../shared_files/tbt_cell_5.mat");
z_binned_cell = importdata("z_binned_cell_ball_norm.mat");
error_cell = importdata("error_cell.mat");
x_binned_cell = importdata("x_binned_cell");
z_correction_corrs_cell = importdata("z_correction_corrs_cell_norm_thresh1.mat");
z_correction_corrs_p_cell = importdata("z_correction_corrs_cell_p_norm_thresh1.mat");

nbins = 50;
ex_md = [1,3];
error_thresh = 1;

% plot_example_neuron_correction_correlation(z_binned_cell,x_binned_cell,error_cell,tbt_cell,z_error_corrs_cell,ex_md);

% [error_cell,x_binned_cell,ro_all,sign_matching] = run_many_check_BMI_correct_errors(virmen_cell,tbt_cell,nbins,false);

[n_ind] = plot_neuron_correction_corr_distribution(z_binned_cell,x_binned_cell,error_cell,tbt_cell,z_correction_corrs_cell,ex_md,error_thresh);

ex_corr = z_correction_corrs_cell{ex_md(1),ex_md(2)}(n_ind);
ex_p = z_correction_corrs_p_cell{ex_md(1),ex_md(2)}(n_ind);

[stats_xc] = plot_neuron_correction_correlations(z_binned_cell,error_cell,tbt_cell,x_binned_cell,error_thresh);

all_res_cell = importdata("all_res_cell_cross_pixels.mat");

% [stats_tests,h_boots] = plot_many_decoding_check(all_res_cell);
[h_boots] = plot_many_corrective_decoding_check(all_res_cell);

Wout_online_cell = importdata("Wout_online_cell.mat");
bmi_weights_cell = importdata("bmi_weights_cell_pixels.mat");

[dot_prods,dot_prods_yaw,h_boots_w] = plot_va_yaw_weights_compare(Wout_online_cell,bmi_weights_cell);