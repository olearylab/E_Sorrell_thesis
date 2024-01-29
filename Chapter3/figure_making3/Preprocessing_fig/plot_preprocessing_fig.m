%% Plot preprocessing figure

% addpath("../shared_functions")

% Load in data
zdata = importdata("DW48_20191119_final.mat");
virmen_data = importdata("trimmed_virmen_DW48_20191119.mat");
results_struct_raw = importdata("DW48_20191119_raw_results_struct_60.mat");
results_struct_dff = importdata("DW48_20191119_dff_results_struct_60.mat");

% Plot example baseline shift in neural activity, shift in decoder output,
% and correction by DF/F filtering.

% Only plotting change in decoder output for now. See if Tim suggests
% adding the neural activity.

plot_DFF_decoder_example(results_struct_raw,results_struct_dff);

% Plot quantitative results
mice = ["DW81";"DW83";"DW86";"DW113";"DW129"];
data_days{1} = ["20200902";"20200907";"20200924"];
data_days{2} = ["20200902";"20200914";"20200922"];
data_days{3} = ["20200902";"20200907";"20200914"];
data_days{4} = ["20210623";"20210628";"20210630";"20210701"];
data_days{5} = ["20211111"];

reps = 5;
training_trials = 80;

[h_boot] = plot_preprocessing_quant(mice,data_days,reps,training_trials);