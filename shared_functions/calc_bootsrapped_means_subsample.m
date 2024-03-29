function [bootstrap_means] = calc_bootsrapped_means_subsample(z_binned,tbt_details,boot_samps,types_vec,sub_sample)
% 07/03/2022

% Function for calculating bootstrapped means of binned neural data.

% I'm expecting this to be very slow. Lots of loops. Maybe need to figure
% out how to speed up.

% resample each trials for each trial type. Whole trials, same for all
% neurons sampled, rather than each bin and/or each neuron independently.

% Do I need anything other than the mean of each resample? e.g. variance?

% 13/01/2023 edit: if not subsampling, bootstrap resample for all trials
% still used smallest number of trials out of all types.

num_types = length(types_vec); % NL, NR, BL, BR. Could change
bootstrap_means = zeros(num_types,boot_samps,size(z_binned,2),size(z_binned,3));

if sub_sample
    [kept_trials] = subsample_trials(tbt_details,types_vec);
    num_samples = size(kept_trials,2);
else
    num_trials = zeros(length(types_vec),1);
    for n = 1:num_types
        cur_trials = find(ismember(tbt_details(3,:),types_vec(n)));
        num_trials(n) = length(cur_trials);
    end
    num_samples = min(num_trials);
end

for n = 1:num_types
    cur_trials = find(ismember(tbt_details(3,:),types_vec(n)));
    if sub_sample
        cur_trials = kept_trials(n,:);
    end
    % num_samples = length(cur_trials);
    for b = 1:boot_samps
        cur_samples = datasample(cur_trials,num_samples);
        cur_resamp = z_binned(cur_samples,:,:);
        bootstrap_means(n,b,:,:) = mean(cur_resamp,'omitnan');
    end
end