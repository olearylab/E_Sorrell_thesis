function [all_means,all_stds] = decoder_dimension_check_n(Wout,Wout_off,zdata,virmen_data,tbt_details,model_params)
% 27/07/2023

% Pass neural activity through decoder weights, as well as through absolute
% decoder weights
% Could maybe pass through only weights in neurons?
% Takes too long to run for all weights, so rerunning for weights in
% neurons

types_vec = [1,4,7,10];

model_params.reg_images = false;
model_params.spatial = false;

[zdatafilt, train_mean, Rfixed] = preprocess_images_05_2021(zdata,model_params,[],[],[],false);

ITI = virmen_data(8,:);
cleaned_valid = clean_valid_data(ITI);

virmen_data = virmen_data(:,cleaned_valid);

trial_num = virmen_data(12,:);

zdatafilt = zdatafilt(cleaned_valid,:);

% These weights are from pixels in neuron.

% Don't include the bias term
on_output = zdatafilt*Wout;
off_output = zdatafilt*Wout_off;

on_output_abs = zdatafilt*abs(Wout);
off_output_abs = zdatafilt*abs(Wout_off);

% Maybe should convert output of yaw decoder to actual yaw? But not really
% sure?
% Should probably subtract off the mean value? Because deviations from the
% mean are what we are after. - have excluded the bias term. This should be
% fine.

ball_t = find(ismember(tbt_details(3,:),types_vec([1,2])));
ball_trials = ismember(trial_num,ball_t);

bmi_t = find(ismember(tbt_details(3,:),types_vec([3,4])));
bmi_trials = ismember(trial_num,bmi_t);

on_ball = mean(abs(on_output(ball_trials)));
on_bmi = mean(abs(on_output(bmi_trials)));
off_ball = mean(abs(off_output(ball_trials)));
off_bmi = mean(abs(off_output(bmi_trials)));

on_ball_s = std(abs(on_output(ball_trials)));
on_bmi_s = std(abs(on_output(bmi_trials)));
off_ball_s = std(abs(off_output(ball_trials)));
off_bmi_s = std(abs(off_output(bmi_trials)));

% on_ball_abs = mean(abs(on_output_abs(ball_trials)));
% on_bmi_abs = mean(abs(on_output_abs(bmi_trials)));
% off_ball_abs = mean(abs(off_output_abs(ball_trials)));
% off_bmi_abs = mean(abs(off_output_abs(bmi_trials)));

on_ball_abs = mean((on_output_abs(ball_trials)));
on_bmi_abs = mean((on_output_abs(bmi_trials)));
off_ball_abs = mean((off_output_abs(ball_trials)));
off_bmi_abs = mean((off_output_abs(bmi_trials)));

on_ball_abs_s = std((on_output_abs(ball_trials)));
on_bmi_abs_s = std((on_output_abs(bmi_trials)));
off_ball_abs_s = std((off_output_abs(ball_trials)));
off_bmi_abs_s = std((off_output_abs(bmi_trials)));

all_means = nan.*ones(4,2);
all_means(1,1) = on_ball;
all_means(1,2) = on_bmi;
all_means(2,1) = off_ball;
all_means(2,2) = off_bmi;
all_means(3,1) = on_ball_abs;
all_means(3,2) = on_bmi_abs;
all_means(4,1) = off_ball_abs;
all_means(4,2) = off_bmi_abs;

all_stds = nan.*ones(4,2);
all_stds(1,1) = on_ball_s;
all_stds(1,2) = on_bmi_s;
all_stds(2,1) = off_ball_s;
all_stds(2,2) = off_bmi_s;
all_stds(3,1) = on_ball_abs_s;
all_stds(3,2) = on_bmi_abs_s;
all_stds(4,1) = off_ball_abs_s;
all_stds(4,2) = off_bmi_abs_s;

% Probably also want to get binned differences...? Come back to this?