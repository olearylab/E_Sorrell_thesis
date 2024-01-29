function [trials_mat,trial_ends,ends_match,trial_corrs] = end_yaw_integration(virmen_data, tbt_details, start_dist, yaw_offset, plot_res)
% 18/11/2021
% function for integrating yaw over end of trial (maybe vary total distance
% over which to integrate). Integrate from current view angle at start of
% integration.

% Could potentially also recreate full trajectory from a point using both
% pitch and yaw.

% There is still an issue in this that the yaw is responding to whatever
% the current view angle is, and I'm not sure there's a way to keep taking
% this into account?

% Only need to worry about BMI trials. Do both correct and incorrect. Want
% to directly compare true trajectories to reconstructed trajectories.
dt = 1/30;
% Clean data
ITI = virmen_data(8,:);
cleaned_valid = clean_valid_data(ITI);
% zdim = size(zdata,2);
trial_num_full = virmen_data(12,:);
% make sure ITI assigned to previous trial
trial_num_full = reassign_ITI_trial_num(cleaned_valid,trial_num_full);
trial_num_clean = trial_num_full(cleaned_valid);
num_trials = max(trial_num_clean);
virmen_data_clean = virmen_data(:,cleaned_valid);

trial_ends = zeros(num_trials,1);
trial_ends_true = zeros(num_trials,1);
trial_corrs = zeros(num_trials,1);
% ends_match = zeros(num_trials,1);

tbt_check = zeros(num_trials,1);
tbt_check(ismember(tbt_details(3,:),[1,5,7,11])) = 1;
tbt_check(ismember(tbt_details(3,:),[2,4,8,10])) = -1;

types_vec = [7,8,10,11];
trials_mat = cell(4,1);
for i = 1:4
    cur_trials = find(tbt_details(3,:) == types_vec(i));
    cur_traj = cell(length(cur_trials),1);
    for j = 1:length(cur_trials)
        cur_trial = virmen_data_clean(:,trial_num_clean==cur_trials(j));
        start_index = find(cur_trial(6,:)>=start_dist,1);
        cur_va = cur_trial(7,start_index:end);
        cur_yaw = cur_trial(15,start_index:end);
        test_valid = ones(length(cur_yaw),1);
        [va] = calc_va_alt_05_2021(cur_yaw,cur_va,yaw_offset,test_valid);
        cur_traj{j} = [cur_va;va];
        trial_ends(cur_trials(j)) = sign(va(end));
        trial_ends_true(cur_trials(j)) = sign(cur_va(end));
        
        % Calculate correlation between integrated VA and decoded VA.
        % Maybe use different metric
        cur_corr = corrcoef(va,cur_va,'rows','pairwise');
        trial_corrs(cur_trials(j)) = cur_corr(1,2);
        
    end
    trials_mat{i} = cur_traj;
end

%ends_match = trial_ends.*tbt_check';
% Currently 1 if matching, -1 if not, 0 if not calculated. Should maybe
% make 0 if not matching.
ends_match = trial_ends.*trial_ends_true;
% ends_match(ends_match~=1) = 0;

if plot_res
    % Visualise trajectories
    titles = ["BMI Left Correct"; "BMI Left Incorrect"; "BMI Right Correct"; "BMI Right Incorrect"]; 
    colour_vec = ["b";"r";"b";"r"];
    figure
    % keep track for bar plot
    bar_vec = zeros(4,1);
    trial_totals = zeros(4,1);
    for i = 1:4
        subplot(2,2,i)
        yline(0,'LineWidth',2)
        title(titles(i))
        hold on
        cur_traj = trials_mat{i};
        trial_totals(i) = length(cur_traj);
        for j = 1:length(cur_traj)
            cur_trials = cur_traj{j};
            plot((1:size(cur_trials,2))*dt,cur_trials(1,:),'LineWidth',2,'Color','k')
            if (sign(cur_trials(1,end))==sign(cur_trials(2,end)))
                plot((1:size(cur_trials,2))*dt,cur_trials(2,:),'LineWidth',2,'Color','b')
                bar_vec(i) = bar_vec(i) + 1;
            else
                plot((1:size(cur_trials,2))*dt,cur_trials(2,:),'LineWidth',2,'Color','r')
            end
            % bar_vec(i) = bar_vec(i) + (sign(cur_trials(1,end))==sign(cur_trials(2,end)));
        end

    end

    % Plot summary of results
    figure
    bar(bar_vec./trial_totals)
end
