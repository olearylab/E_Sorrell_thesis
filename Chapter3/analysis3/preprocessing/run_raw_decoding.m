%% Run without preprocessing

initialise_params;
% Can make these all true for with preprocessing
model_params.reg_images = false;
model_params.spatial = false;
model_params.dff = false;

mice = ["DW81";"DW83";"DW86"];
data_days{1} = ["20200902";"20200907";"20200924"];
data_days{2} = ["20200902";"20200914";"20200922"];
data_days{3} = ["20200902";"20200907";"20200914"];

% Did I change the dff_offset for the offline analysis? I don't think so

run_save_offline(mice, data_days, model_params, 5, 80, true);

mice = ["DW113";"DW129"];
clear data_days
data_days{1} = ["20210623";"20210628";"20210630";"20210701"];
data_days{2} = ["20211111"];

run_save_offline_2021_mice(mice, data_days, model_params, 5, 80, true);