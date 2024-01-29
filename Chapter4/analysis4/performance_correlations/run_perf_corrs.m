%% Run analysis for performance correlations

% Just need to get hp_results_all, performance_vec and lp_RMSE_all

% This was orginally run on more than the actual data used for analysis.
% This is trimmed in the plotting code.

% load in the full cells

% What are var_samps????
var_samps = 31;
hp_filter_cut = 1;
yaw_offsets = [1.495.*ones(23,1);1.4844.*ones(17,1)];
mice_vec = [81;81;81;81;81;81;81;81;83;83;83;83;83;83;83;86;86;86;86;86;86;86;86;113;113;113;113;113;113;113;111;111;111;116;116;116;129;129;129;129];
[hp_results_all,yaw_results_all,performance_vec] = plot_variability_performance_correlations(virmen_cell_full,tbt_cell_full,summary_cell_full,mice_vec,yaw_offsets,hp_filter_cut,var_samps);

% I assume this is what I did but I should rerun.
[RMSE_all,R_square_all] = calc_many_online_decoder_errors(virmen_cell_full,tbt_cell_full);