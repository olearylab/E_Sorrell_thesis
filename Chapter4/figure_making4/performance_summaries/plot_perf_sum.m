%% FIGURE 2
addpath("../../../shared_functions")
% prepare workspace for plotting
% macbook
% res_path = "/Users/ethansorrell/Dropbox (Cambridge University)/Boston_PPC_Data_05_2021/";
% big mac
% res_path = "/Users/ethan/Dropbox (Cambridge University)/Boston_PPC_Data_05_2021/"; 
virmen_cell = importdata("../../../shared_files/virmen_cell_5.mat");
tbt_cell = importdata("../../../shared_files/tbt_cell_5.mat");

trial_types_summary_cell = importdata("trial_types_summary_cell_w_train_5.mat");
days_cell{1} = [0,1,2,3,4];
days_cell{2} = [0,1,2,3,4];
days_cell{3} = [0,1,2,3,4];
days_cell{4} = [0,1,3,4];
days_cell{5} = [0,1,2,3,4];

virmen_cell{4,2} = [];

% And plot

[sig_results,p,summary_results] = plot_BMI_performance_all_mice_w_train_06122022(trial_types_summary_cell,days_cell);

summary_cell = importdata("../../../shared_files/t_types_summary_cell_5.mat");
plot_perf_across_days(summary_cell);