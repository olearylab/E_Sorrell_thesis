%%% Plot weights analysis 
%%
% 
res1 = importdata("DW81_0902_results_struct_11042023.mat");
% results_struct = importdata("res_struct_DW113_20210701_rep1_80.mat");
% ex_weights = results_struct.Wout;
ex_weights = res1.Wout;
stat_cell_CNN_offline = importdata("stat_cell_CNN_offline.mat");
% ex_stat = stat_cell_CNN_offline{13};
ex_stat = stat_cell_CNN_offline{1};

% ex_im = importdata("example_mean_training_im.mat");
ex_im = importdata("DW81_0902_mean_train_im_trans.mat");
trans_weights = true;
plot_example_weights(ex_weights,ex_stat,ex_im,trans_weights);


Wout = res1.Wout(:,1);
full_all_res = importdata("full_all_res_DW81_0902_11042023_ypos.mat");
[thresh_num,thresh_num_RMSE, prop_90] = plot_decoder_dimensionality_and_n_weights_11042023(full_all_res,Wout,0.9,stat_cell_CNN_offline{1});

mice = ["DW81";"DW83";"DW86";"DW113";"DW129"];
data_days{1} = ["20200902";"20200907";"20200924"];
data_days{2} = ["20200902";"20200914";"20200922"];
data_days{3} = ["20200902";"20200907";"20200914"];
data_days{4} = ["20210623";"20210628";"20210630";"20210701"];
data_days{5} = ["20211111"];

rep = 1;
training_trials = 80;
[all_corrs,corr_val,p_val,h_boots] = compare_pix_neu_many(mice,data_days,rep,training_trials,stat_cell_CNN_offline);

