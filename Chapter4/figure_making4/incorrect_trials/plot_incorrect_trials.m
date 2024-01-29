%% FIGURE 7
% addpath("../shared_functions")
% prepare workspace for plotting figure 7
% macbook
% res_path = "/Users/ethansorrell/Dropbox (Cambridge University)/Boston_PPC_Data_05_2021/";
% big mac
% res_path = "/Users/ethan/Dropbox (Cambridge University)/Boston_PPC_Data_05_2021/"; 
virmen_cell = importdata("../../../shared_files/virmen_cell_5.mat");
tbt_cell = importdata("../../../shared_files/tbt_cell_5.mat");

virmen_cell{4,2} = [];

virmen_data = virmen_cell{4,1};
t_num = 23;
start_dist = 295;
yaw_offset = 1.4844;
plot_example_yaw_integration_09012023(virmen_data,t_num,start_dist,yaw_offset);

% results_struct_cell = importdata("results_struct_cell_5.mat");
results_struct_cell = importdata("results_struct_cell_pixels.mat");
% ind = 3;
ind = 1;
[h_boots] = position_decoding_summary_plot_09012023(results_struct_cell, virmen_cell, tbt_cell, ind);

virmen_data = virmen_cell{5,3};
tbt_details = tbt_cell{5,3};
plot_res = true;
[trials_mat,trial_ends,ends_match,trial_corrs] = end_yaw_integration_09012023(virmen_data, tbt_details, start_dist, yaw_offset, plot_res);

offset_vec = [1.4961,1.4961,1.4961,1.4844,1.4844];
% end_traces_match_scores_09012023(virmen_cell,tbt_cell,start_dist,offset_vec);
% [h_boots] = end_traces_match_scores_12062023(virmen_cell,tbt_cell,start_dist,offset_vec);
[h_boots] = end_traces_match_scores_11122023(virmen_cell,tbt_cell,start_dist,offset_vec);

nbins =50;
classifiers_res_cell = importdata("classifiers_res_cell_svm_5.mat");
neural_classifiers_summary_analysis_09012023(classifiers_res_cell,virmen_cell,tbt_cell,nbins);

[conf_mat, full_compare, full_compare_correct] = compare_behaviour_neural_classification_09012023(classifiers_res_cell,results_struct_cell,virmen_cell,tbt_cell,nbins,start_dist,offset_vec);