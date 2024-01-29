function [all_processed,all_vars,final_results,final_yaw_results] = trial_averaged_variance_measures_full(virmen_data,tbt_details,var_samps,hp_cut,mean_samps,yaw_offset)
% 16/05/22

% Compare measures of local view angle variance. Calculate for individual
% trials then average across trials. Hopefully this removes issue of
% discontinuity across trials. Maybe also need to remove x samples from
% beginning and end of each trial aswell?

% Look at moving variance, variance of trial - moving average, and variance
% of high passed version of trial.

% (signal - moving average is surely just like signal - low pass filter? In
% which case might as well be using high pass filtered version)
% Realistically want to identify frequency present in decoded output and
% not in true view angle.

% Edited from previous version to include open-loop va decoding measures as
% well as yaw moving average measure (maybe pick a better yaw measure).


fs = 30;

%%
circum = 64;
V = 0.32;
alpha = -50/75*circum/V;
beta = 0.05*circum/V/2.5;
dt = 1/30;

va_vel = yaw_offset + (virmen_data(17,2:end)-virmen_data(17,1:end-1))./(dt*-beta);
va_vel = [va_vel,0];
virmen_data = [virmen_data;va_vel];

%%
ITI = virmen_data(8,:);
cleaned_valid = clean_valid_data(ITI);
virmen_data = virmen_data(:,cleaned_valid);
% Don't want to wrap view angle as this creates discontinuities. 

%% 
control_trials = find(ismember(tbt_details(3,:),[1,4]));
bmi_trials = find(ismember(tbt_details(3,:),[7,10]));
correct_trials = find(ismember(tbt_details(3,:),[1,4,7,10]));

% moving variance, high pass filter, signal - moving average.
% view angle, decoded view angle
all_processed = zeros(3,2,size(virmen_data,2));
all_vars = zeros(3,2,size(tbt_details,2));
% mov_vars = zeros();
% hp_var = zeros();
% minus_av_var = zeros();

% Just moving average for yaw. Maybe change later
yaw_processed = zeros(2,size(virmen_data,2));
yaw_vars = zeros(2,size(tbt_details,2));

trial_num = virmen_data(12,:);

for n = 1:size(tbt_details,2)
    cur_trial = virmen_data(7,trial_num==n);
    all_processed(1,1,trial_num==n) = movvar(cur_trial,var_samps);
    all_processed(2,1,trial_num==n) = highpass(cur_trial,hp_cut,fs);
    all_processed(3,1,trial_num==n) = cur_trial - movmean(cur_trial,mean_samps);
    
    all_vars(1,1,n) = mean(all_processed(1,1,trial_num==n));
    all_vars(2,1,n) = var(all_processed(2,1,trial_num==n));
    all_vars(3,1,n) = var(all_processed(3,1,trial_num==n));
    
    cur_trial = virmen_data(17,trial_num==n);
    all_processed(1,2,trial_num==n) = movvar(cur_trial,var_samps);
    all_processed(2,2,trial_num==n) = highpass(cur_trial,hp_cut,fs);
    all_processed(3,2,trial_num==n) = cur_trial - movmean(cur_trial,mean_samps);
   
    all_vars(1,2,n) = mean(all_processed(1,2,trial_num==n));
    all_vars(2,2,n) = var(all_processed(2,2,trial_num==n));
    all_vars(3,2,n) = var(all_processed(3,2,trial_num==n));
    
    cur_yaw_trial = virmen_data(15,trial_num==n);
   	yaw_processed(1,trial_num==n) = movvar(cur_yaw_trial,var_samps);
    
    yaw_vars(1,n) = mean(yaw_processed(1,trial_num==n));
    
    cur_yaw_trial = virmen_data(end,trial_num==n);
   	yaw_processed(2,trial_num==n) = movvar(cur_yaw_trial,var_samps);
    
    yaw_vars(2,n) = mean(yaw_processed(2,trial_num==n));
    
end

% process type x control/bmi x observed/decoded
final_results = zeros(3,2,2);
for i = 1:2
    final_results(1,1,i) = mean(all_vars(1,i,control_trials));
    final_results(1,2,i) = mean(all_vars(1,i,bmi_trials));
    final_results(2,1,i) = mean(all_vars(2,i,control_trials));
    final_results(2,2,i) = mean(all_vars(2,i,bmi_trials));
    final_results(3,1,i) = mean(all_vars(3,i,control_trials));
    final_results(3,2,i) = mean(all_vars(3,i,bmi_trials));
end

% control/bmi x observed/decoded
final_yaw_results = zeros(2,2);
for i = 1:2
    final_yaw_results(1,i) = mean(yaw_vars(i,control_trials));
    final_yaw_results(2,i) = mean(yaw_vars(i,bmi_trials));
end
