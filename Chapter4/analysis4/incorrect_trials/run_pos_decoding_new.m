%% Run offline decoding of online data
initialise_params;
% Also need to create the weights and train mean cells. Just need to load
% in the training data for each mouse and run
[Wout, xpred, origztrainfilt, train_mean, Rfixed, train_valid, model_params] = DW_online_training_05_2021(ztrain, virmen_data, model_params);

mice = ["DW81";"DW83";"DW86"];
data_days{1} = ["20200902";"20200907";"20200924"];
data_days{2} = ["20200902";"20200914";"20200922"];
data_days{3} = ["20200902";"20200907";"20200914"];

% Did I change the dff_offset for the offline analysis? I don't think so

run_save_online(mice, data_days, model_params, Wout_train_cell, train_mean_cell, Rfixed_cell, true);

mice = ["DW113";"DW129"];
clear data_days
data_days{1} = ["20210623";"20210628";"20210630";"20210701"];
data_days{2} = ["20211111"];

run_save_online_2021_mice(mice, data_days, model_params, Wout_train_cell, train_mean_cell, Rfixed_cell, true);