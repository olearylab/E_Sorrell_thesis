function [all_trials_res,trials_summary_means,trials_summary_stds] = trial_by_trial_accuracy(results_struct, xfull, tbt_details)

% assess accuracy of decoding on a trial by trial basis.
% want accuracy for every individual trial.
% also return statistics for each trial type.

% Possibly need to take into account trials where mouse starts disengaged
% but still finishes the trial (i.e. the turned greater than 2 pi trials).
test_valid = results_struct.test_valid;
xtest = results_struct.xtest(test_valid,:);
xprediction = results_struct.xprediction(test_valid,:);
trial_num = xfull(12,test_valid);

num_x = size(xtest,2);

all_trials_res = zeros(size(tbt_details,2),3,num_x);

for k = 1:size(tbt_details,2)

    RMSE = sqrt(mean((xtest(trial_num==k,:) - xprediction(trial_num==k,:)).^2,1));
    r_p = zeros(1,num_x);
    R_square = zeros(1,num_x);
    for i = 1:num_x
        cur_cor = corrcoef(xtest(trial_num==k,i),xprediction(trial_num==k,i));
        r_p(i) = cur_cor(1,2);
        R_square(i) = 1 - sum((xtest(trial_num==k,i) - xprediction(trial_num==k,i)).^2)/sum((xtest(trial_num==k,i) - mean(xtest(trial_num==k,i))).^2);
    end
    
    all_trials_res(k,:,:) = [RMSE;R_square;r_p];
    
end

% calculate mean and std of accuracy for each trial type

trials_summary_means = zeros(12,3,num_x);
trials_summary_stds = zeros(12,3,num_x);
for i = 1:12
    % works along first dimension by default
    if sum(tbt_details(3,:)==i) == 1
        trials_summary_means(i,:,:) = all_trials_res(tbt_details(3,:)==i,:,:);
        trials_summary_stds(i,:,:) = 0;
    else
        trials_summary_means(i,:,:) = mean(all_trials_res(tbt_details(3,:)==i,:,:));
        trials_summary_stds(i,:,:) = std(all_trials_res(tbt_details(3,:)==i,:,:));
    end
end
    