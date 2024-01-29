% Plot decoder comparisons figure

addpath("../shared_files")

mice = ["DW81";"DW83";"DW86";"DW113";"DW129"];
data_days{1} = ["20200902";"20200907";"20200924"];
data_days{2} = ["20200902";"20200914";"20200922"];
data_days{3} = ["20200902";"20200907";"20200914"];
data_days{4} = ["20210623";"20210628";"20210630";"20210701"];
data_days{5} = ["20211111"];

example_m_d = "DW113_20210701";
reps = 5;
training_trials = 80;

[h_boots] = plot_offline_decoder_comparisons(mice,data_days,example_m_d,reps,training_trials);