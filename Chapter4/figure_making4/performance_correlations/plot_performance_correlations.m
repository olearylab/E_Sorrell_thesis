%% FIGURE 4
% addpath("../shared_functions")
% Prepare workspace for Figure 4 plotting
% macbook
% res_path = "/Users/ethansorrell/Dropbox (Cambridge University)/Boston_PPC_Data_05_2021/";
% big mac
% res_path = "/Users/ethan/Dropbox (Cambridge University)/Boston_PPC_Data_05_2021/"; 
virmen_cell = importdata("../../../shared_files/virmen_cell_5.mat");
tbt_cell = importdata("../../../shared_files/tbt_cell_5.mat");
summary_cell = importdata("../../../shared_files/t_types_summary_cell_5.mat");
hp_results_all = importdata("hp_results_all_1hz.mat");
RMSE_all = importdata("lp_RMSE_all.mat");
plot_inds = importdata("plot_inds.mat");
performance_vec = importdata("performance_vec.mat");
mice_inds = importdata("mice_inds.mat");

nbins = 50;
linearise_x = true;
ex_offsets = [1.4953,1.4953,1.4953,1.4804,1.4804];
ex_mice = [1,4];

ex_mds = [1,1;4,1];
ex_trials = [0,25;0,21];
hp_cut = 1;

% And plot

[RMSE_temp,corr_vals,p_vals] = plot_biases_and_correlations_06122022(virmen_cell,tbt_cell,summary_cell,nbins,linearise_x,ex_offsets,ex_mice);

[corr_vals2,p_vals2] = plot_va_hp_variance_analysis_combined_06122022(ex_mds,ex_trials,virmen_cell, hp_results_all, performance_vec, plot_inds, mice_inds, hp_cut, RMSE_all);
