function [error_mat,x_binned] = check_error_correction_07032023(virmen_data,tbt_details,nbins,offsets,plot_res)
% 07/03/2023

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

mean_binned = nan.*ones(size(x_binned,2),size(x_binned,3),4);
for i = 1:4
    mean_binned(:,:,i) = squeeze(mean(x_binned(tbt_details(3,:)==types_vec(i),:,:),1,'omitnan'));
end

% Calculate error signal

error_mat = nan.*ones(size(x_binned,1),size(x_binned,2));
% Edited to also calculate for ball trials
mean_ind = [1,2,1,2];
for i = 1:4
    cur_trials = find(tbt_details(3,:)==types_vec(i));
    for n = 1:length(cur_trials)
        error_mat(cur_trials(n),:) = squeeze(mean_binned(:,2,mean_ind(i))) - squeeze(x_binned(cur_trials(n),:,2))';
    end
end

if plot_res
    figure

    for i = 1:2
        subplot(1,2,i)
        cur_trials = find(tbt_details(3,:)==types_vec(2+i));
        cur_errors = error_mat(cur_trials,:);
        cur_yaw = squeeze(x_binned(cur_trials,:,3));
        ro = corr(cur_errors(:),cur_yaw(:));
        scatter(cur_errors(:),cur_yaw(:));
        title(ro)
    end
end