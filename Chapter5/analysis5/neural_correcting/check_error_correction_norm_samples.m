function [error_mat,correcting_mat] = check_error_correction_norm_samples(virmen_data,tbt_details,nbins,offsets,mean_binned,std_binned)
% 08/09/2023

% Function for assessing whether physical "yaw" during bmi trials could be
% an error correcting signal?

% Am I sure about the sign of yaw vs view angle? Do I need to change
% anything for this.
types_vec = [1,4,7,10];
inc_vec = [2,5,8,11]; % Incorrect trials

circum = 64;
V = 0.32;
alpha = -50/75*circum/V;
beta = 0.05*circum/V/2.5;
virmen_data(13,:) = alpha*(virmen_data(13,:)-offsets(1));
virmen_data(15,:) = -beta*(virmen_data(15,:)-offsets(2));

linearise_x = true;
color_vec = [[0.4940 0.1840 0.5560];[0.4660 0.6740 0.1880]];
x_vec = [6,7,15,17]; % Changed from [6,7,15];
virmen_data(7,:) = wrapToPi(virmen_data(7,:));
[x_binned,centres] = bin_kin_data(virmen_data,x_vec,linearise_x,nbins);

% mean_binned = nan.*ones(size(x_binned,2),size(x_binned,3),4);
% std_binned = nan.*ones(size(x_binned,2),size(x_binned,3),4);
% for i = 1:4
%     mean_binned(:,:,i) = squeeze(mean(x_binned(tbt_details(3,:)==types_vec(i),:,:),1,'omitnan'));
%     std_binned(:,:,i) = squeeze(std(x_binned(tbt_details(3,:)==types_vec(i),:,:),0,1,'omitnan'));
% end

% Calculate what bin samples are in (using same edges as mean va)
ITI = virmen_data(8,:);
cleaned_valid = clean_valid_data(ITI);
virmen_data(6,:) = virmen_data(6,:) + abs(virmen_data(5,:));
virmen_data_clean = virmen_data(:,cleaned_valid);
[binned,edges] = discretize(virmen_data_clean(6,:),nbins);
[binned,edges] = discretize(virmen_data(6,:),edges);
trial_num = virmen_data(12,:);

% Calculate error signal. Normalise by va std

error_mat = nan.*ones(1,size(virmen_data,2));
% Edited to also calculate for ball trials
mean_ind = [1,2,1,2];
for i = 1:4
    cur_trials = find(tbt_details(3,:)==types_vec(i));
    cur_trial_num = ismember(trial_num,cur_trials);
    for b = 1:nbins
        error_mat(1,(binned==b)&cur_trial_num) = (virmen_data(7,(binned==b)&cur_trial_num) - squeeze(mean_binned(b,2,mean_ind(i))))./squeeze(std_binned(b,2,mean_ind(i)));
    end
end

% +ve error requires -ve view angle velocity to be correcting
correcting_mat = sign(virmen_data(15,:)) ~= sign(error_mat);