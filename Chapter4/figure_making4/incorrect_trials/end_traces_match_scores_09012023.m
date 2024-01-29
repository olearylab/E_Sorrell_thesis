function [] = end_traces_match_scores_09012023(virmen_cell,tbt_cell,start_dist,offset_vec)
% function for running end yaw integration for many days and many mice. 

% Need to come up with a score for each mouse to compare across mice.

% Loop through mice and days, calculating proportion matching end results
% for each type, combining days together. Potentially combining left and
% right together as well.

% Tim also suggested to be careful about normalisation due to different
% lengths of trials from start distance. I guess he is suggesting that over
% time deviations accumulate, thus simple mean is not sufficient. Need to
% think about this

num_mice = size(virmen_cell,1);
num_days = size(virmen_cell,2);

%%
% TODO: Save correlations for each mouse. collect each day and type
% separately and can combine however necessary for plotting. 
corr_results = cell(num_mice,num_days);
corr_results_combine_days = cell(num_mice,1);
ends_match_combine_days = cell(num_mice,1);
ends_match_cell = cell(num_mice,num_days);
full_m_tbt = cell(num_mice,1);
for m = 1:num_mice
    cur_m_tbt = [];
    cur_m_corrs = [];
    cur_m_ends = [];
    for d = 1:num_days
        if ~isempty(virmen_cell{m,d})
            [trials_mat,trial_ends,ends_match,trial_corrs] = end_yaw_integration(virmen_cell{m,d}, tbt_cell{m,d}, start_dist, offset_vec(m), false);
            % trials_mat is 4x1 cell, containing cell with each trial for blc,
            % bli, brc, bri. Each is 2xt true va and recovered va.
            % Would need to change function to get normal trials (which would
            % have quite high correlations).

            % Store correlations
            % trial_corrs is one value for each trial
            corr_results{m,d} = trial_corrs;  
            ends_match_cell{m,d} = ends_match;
            cur_m_corrs = [cur_m_corrs;trial_corrs];
            cur_tbt = tbt_cell{m,d};
            cur_m_tbt = [cur_m_tbt,cur_tbt(3,:)];
            cur_m_ends = [cur_m_ends;ends_match];
        %else
            %corr_results{m,d} = nan;  
            %ends_match_cell{m,d} = nan;
        end
    end
    full_m_tbt{m} = cur_m_tbt;
    corr_results_combine_days{m} = cur_m_corrs;
    ends_match_combine_days{m} = cur_m_ends;
end

%% Visualisation

% plot distribution of scores, where each trial is scored independently
% Plot mice together.
% create long vector and separate into groups.
% figure
% types_vec = [7,8,10,11];
% full_res = cell(4,1);
% full_groups = cell(4,1);
% for m = 1:num_mice
%     cur_m = corr_results_combine_days{m};
%     cur_tbt = full_m_tbt{m};
%     for j = 1:4
%         % blc, bli, brc, bri
%         full_res{j} = [full_res{j};cur_m(cur_tbt==types_vec(j))];
%         full_groups{j} = [full_groups{j};m*ones(sum(cur_tbt==types_vec(j)),1)];
%     end
% end
% 
% titles = ["BMI Left Correct";"BMI Left Incorrect";"BMI Right Correct";"BMI Right Incorrect"];
% for j = 1:4
%     subplot(1,4,j)
%     yline(0,'LineWidth',2);
%     hold on
%     % blc, bli, brc, bri
%     boxplot(full_res{j},full_groups{j});
%     title(titles(j))
%     xlabel("Mouse")
%     ylim([-1,1])
%     xticks(1:num_mice)
%     if j == 1
%         ylabel("Correlations")
%         yticks([-1,-0.5,0,0.5,1])
%     else
%         yticks([])
%     end
% end
%


% Plot mice separated.
% figure
% types_vec = [7,8,10,11];
% full_res = cell(num_mice,1);
% full_groups = cell(num_mice,1);
% for m = 1:num_mice
%     cur_m = corr_results_combine_days{m};
%     cur_tbt = full_m_tbt{m};
%     for j = 1:4
%         % blc, bli, brc, bri
%         full_res{m} = [full_res{m};cur_m(cur_tbt==types_vec(j))];
%         full_groups{m} = [full_groups{m};j*ones(sum(cur_tbt==types_vec(j)),1)];
%     end
% end
% 
% for m = 1:num_mice
%     subplot(1,num_mice,m)
%     yline(0,'LineWidth',2);
%     hold on
%     boxplot(full_res{m},full_groups{m});
%     title("Mouse " + m)
%     xlabel("Trial Type")
%     ylim([-1,1])
%     xticks(1:4)
%     if m == 1
%         ylabel("Correlations")
%         yticks([-1,-0.5,0,0.5,1])
%     else
%         yticks([])
%     end
%     
% end

% Combine left and right
figure
types_mat = [7,10;8,11];
full_res = cell(2,1);
full_groups = cell(2,1);
for m = 1:num_mice
    cur_m = corr_results_combine_days{m};
    cur_tbt = full_m_tbt{m};
    for j = 1:2
        % blc, bli, brc, bri
        full_res{j} = [full_res{j};cur_m(ismember(cur_tbt,types_mat(j,:)))];
        full_groups{j} = [full_groups{j};m*ones(sum(ismember(cur_tbt,types_mat(j,:))),1)];
    end
end
titles_2 = ["BMI Correct","BMI Incorrect"];
for j = 1:2
    subplot(1,2,j)
    yline(0,'LineWidth',2);
    hold on
    % blc, bli, brc, bri
    boxplot(full_res{j},full_groups{j});
    title(titles_2(j));
    xlabel("Mouse")
    ylim([-1,1])
    xticks(1:num_mice)
    if j == 1
        ylabel("Correlations")
        yticks([-1,-0.5,0,0.5,1])
    else
        yticks([])
    end
end


% Look for trend over time?
% Still unsure if mouse together or mouse separate is better
% 
% figure
% types_vec = [7,8,10,11];
% full_means = zeros(num_mice,num_days,4);
% full_stds = zeros(num_mice,num_days,4);
% for m = 1:num_mice
%     for d = 1:num_days
%         cur_m = corr_results{m,d};
%         cur_tbt = tbt_cell{m,d};
%         for j = 1:4
%             % blc, bli, brc, bri
%             if ~isempty(cur_m)
%                 full_means(m,d,j) = mean(cur_m(cur_tbt(3,:)==types_vec(j)),'omitnan');
%                 full_stds(m,d,j) = std(cur_m(cur_tbt(3,:)==types_vec(j)),'omitnan');
%             else
%                 full_means(m,d,j) = nan;
%                 full_stds(m,d,j) = nan;
%             end
%         end
%     end
% end
% 
% for j = 1:4
%     subplot(1,4,j)
%     hold on
%     % blc, bli, brc, bri
%     errorbar(squeeze(full_means(:,:,j))',squeeze(full_stds(:,:,j))','LineWidth',2)
%     title(titles(j))
%     xlabel("Day")
%     ylim([-1,1])
%     xticks(1:num_days)
%     if j == 1
%         ylabel("Mean Correlations")
%         yticks([-1,-0.5,0,0.5,1])
%     else
%         yticks([])
%     end
% end


% Compare matching sign of decoder and integrated yaw at end.
% Potentially add numbers of trials to plots.
% Potentially summarise both directions together.
% Mice Together
% figure
% types_vec = [7,8,10,11];
% bar_mat = zeros(num_mice,4);
% for m = 1:num_mice
%     cur_m = ends_match_combine_days{m};
%     cur_tbt = full_m_tbt{m};
%     for j = 1:4
%         % blc, bli, brc, bri
%         bar_mat(m,j) = sum(cur_m(cur_tbt==types_vec(j))==1)/sum(cur_tbt==types_vec(j));
%     end
% end
% 
% titles = ["BMI Left Correct";"BMI Left Incorrect";"BMI Right Correct";"BMI Right Incorrect"];
% for j = 1:4
%     subplot(2,2,j)
%     % yline(0,'LineWidth',2);
%     hold on
%     % blc, bli, brc, bri
%     bar(bar_mat(:,j),'w','LineWidth',2)
%     title(titles(j))
%     if j >2
%         xlabel("Mouse")
%     end
%     ylim([0,1])
%     xticks(1:num_mice)
%     if (j == 1) || (j == 3)
%         ylabel("Proportion Matching")
%         yticks([0,0.5,1])
%     else
%         yticks([])
%     end
% end

% Mice Separate
% figure
% for m = 1:num_mice
%     subplot(1,num_mice,m)
%     % yline(0,'LineWidth',2);
%     hold on
%     % blc, bli, brc, bri
%     bar(bar_mat(m,:),'w','LineWidth',2)
%     title("Mouse " + m)
%     xlabel("Trial Type")
%     ylim([0,1])
%     xticks(1:4)
%     if (m == 1) 
%         ylabel("Proportion Matching")
%         yticks([0,0.5,1])
%     else
%         yticks([])
%     end
% end

% Combine left and right

% Mice Together
figure
types_vec = [7,10;8,11];
bar_mat = zeros(num_mice,4);
for m = 1:num_mice
    cur_m = ends_match_combine_days{m};
    cur_tbt = full_m_tbt{m};
    for j = 1:2
        % blc, bli, brc, bri
        bar_mat(m,j) = sum(cur_m(ismember(cur_tbt,types_vec(j,:)))==1)/sum(ismember(cur_tbt,types_vec(j,:)));
    end
end

titles = ["BMI Correct";"BMI Incorrect"];
for j = 1:2
    subplot(1,2,j)
    % yline(0,'LineWidth',2);
    hold on
    % blc, bli, brc, bri
    bar(bar_mat(:,j),'w','LineWidth',2)
    title(titles(j))
    xlabel("Mouse")
    ylim([0,1])
    xticks(1:num_mice)
    if (j == 1)
        ylabel("Proportion Matching")
        yticks([0,0.5,1])
    else
        yticks([])
    end
end