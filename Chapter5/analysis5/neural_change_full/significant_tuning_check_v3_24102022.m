function [bootstrap_means,sig_results,sig_check,centres] = significant_tuning_check_v3_24102022(zdata, virmen_data, tbt_details, num_shuffles, linearise_x, nbins, boot_samps, shuff_limit, xnum)
% 24/10/2022

% Function for determining cells with significant tuning.
% Using circular shuffle check.
% Maybe use SNR comparison between unshuffled and shuffled as significance
% check.

% Dan says significance check is whether max value is higher than 99% of
% shuffles max (potentially more fine check than that).

% In this case there is no real connection to trial type anymore, so bin
% all shuffles as single curve, then no bootstrapping.
% Alternatively still bin into trials, but all one trial type.

% Do I remove ITI before or after shuffling? Seems odd to include invalid
% data, but continuity in time removed if removing invalid first.

% Will likely take a long time.

% 24/10/2022 change: Nan out all BMI trial timepoints

zdim = size(zdata,2);
t_types = [1,4,7,10];


%% remove invalid data
ITI = virmen_data(8,:);
cleaned_valid = clean_valid_data(ITI);

trial_num_full = virmen_data(12,:);
% make sure ITI assigned to previous trial
trial_num_full = reassign_ITI_trial_num(cleaned_valid,trial_num_full);

% Not removing incorrect trials as invalid (this was possibly a bug in
% previous versions). Nan out later if desired.
% correct_trials = find(ismember(tbt_details(3,:),t_types));
% cleaned_valid = cleaned_valid & ismember(trial_num_full,correct_trials);

trial_num_clean = trial_num_full(cleaned_valid);

virmen_data_clean = virmen_data(:,cleaned_valid);
%% Preprocess and bin data

initialise_params;
if ndims(zdata) == 2
    model_params.spatial = false;
    create_reg = false;
else
    create_reg = true;
end
model_params.reg_images = false;
disp("Pre-Processing")
[zdatafilt, train_mean, Rfixed] = preprocess_images_05_2021(zdata,model_params,[],[],[],create_reg);
disp("Pre-Processing Complete")   

% Need to change if I want to zscore properly
zdatafilt_clean = zdatafilt(cleaned_valid,:);

num_samples = size(zdata,1);

%% Nan out bmi trials
zdatafilt(ismember(virmen_data(1,:),[3,4]),:) = nan;

%% Perform shuffles
% Maybe make shuffle range an input
shuff_inds = randi([shuff_limit,(num_samples-shuff_limit)],num_shuffles,1);
% shuff_bootstrap_means = zeros(length(t_types),num_shuffles,nbins,zdim);
% shuff_bootstrap_means_max = zeros(length(t_types),num_shuffles,zdim);

% shuff_z_binned_means_max = zeros(length(t_types),num_shuffles,zdim);
shuff_types = [1,4];
shuff_z_binned_means_max = zeros(length(shuff_types),num_shuffles,zdim);
% Shuffle using full data (not cleaned data).
for i = 1:num_shuffles
    shuff_zdatafilt = circshift(zdatafilt,shuff_inds(i));
    
    % Bin into single curve
    % z_binned = nbins x zdim
    % [z_binned, centres] = bin_neural_data_general_all(virmen_data,shuff_zdatafilt,linearise_x,nbins);
    
    % Bin into trials
    [z_binned, centres] = bin_neural_data_general_any(virmen_data_clean,shuff_zdatafilt(cleaned_valid,:),xnum,linearise_x,nbins);
    % Potentially bootstrap this to get estimate of mean shuffled tuning
    
    %z_binned = squeeze(mean(z_binned,'omitnan'));
    z_binned_means = zeros(length(shuff_types),nbins,zdim);
    for j = 1:length(shuff_types)
        % cur_trials = find(tbt_details(3,:)==t_types(j));
        z_binned_means(j,:,:) = mean(z_binned(tbt_details(3,:)==shuff_types(j),:,:),'omitnan');
    end
    
    
    % Bootstrap
    % bootstrap means are num_types x bootsampes x nbins x zdim
    % [bootstrap_means] = calc_bootsrapped_means(z_binned,tbt_details,boot_samps,t_types);
    
    % Calculate for each type, then maximum over all types
    % cur_mean = squeeze(mean(bootstrap_means,2,'omitnan'));
    % shuff_bootstrap_means_max(:,i,:) = squeeze(max(cur_mean,[],2,'omitnan'));
    
    shuff_z_binned_means_max(:,i,:) = squeeze(max(z_binned_means,[],2,'omitnan'));
    
end
shuff_means_max = squeeze(max(shuff_z_binned_means_max,[],1,'omitnan'));
%% Unshuffled
% z_binned is trials x nbins x zdim
[z_binned, centres] = bin_neural_data_general_any(virmen_data_clean,zdatafilt_clean,xnum,linearise_x,nbins);

% Use bootstrapping to calculate distributions of sample means for each
% trial type in each bin. 
% Number of resamples defined by boot_samps
% bootstrap_means is num_types x boot_samps x bins x neurons

[bootstrap_means] = calc_bootsrapped_means(z_binned,tbt_details,boot_samps,t_types);
    
%% Determine significant tuning
sig_results = zeros(length(t_types),zdim);

% Potentially add confidence interval on maximum bootstrap mean. Check if
% lower 1% rather than mean

boot_means_mean = squeeze(mean(bootstrap_means,2,'omitnan'));
boot_means_mean_max = squeeze(max(boot_means_mean,[],2,'omitnan'));

sig_check = squeeze(prctile(shuff_means_max,99,1));

for i = 1:4
    sig_results(i,:) = boot_means_mean_max(i,:) > sig_check;
end