function [h_boots] = end_traces_match_scores_11122023(virmen_cell,tbt_cell,start_dist,offset_vec)
% 11/12/2023
% function for running end yaw integration for many days and many mice. 

% Need to come up with a score for each mouse to compare across mice.

% Loop through mice and days, calculating proportion matching end results
% for each type, combining days together. Potentially combining left and
% right together as well.


% Edited to lump mice together
marker_size = 10;

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

mean_corr_results = nan.*ones(num_mice,num_days,2);
mean_ends_match = nan.*ones(num_mice,num_days,2);
types_mat = [7,10;8,11];

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

            mean_corr_results(m,d,1) = mean(trial_corrs(ismember(cur_tbt(3,:),types_mat(1,:))));
            mean_corr_results(m,d,2) = mean(trial_corrs(ismember(cur_tbt(3,:),types_mat(2,:))));

            mean_ends_match(m,d,1) = sum(ends_match(ismember(cur_tbt(3,:),types_mat(1,:)))==1)/sum(ismember(cur_tbt(3,:),types_mat(1,:)));
            mean_ends_match(m,d,2) = sum(ends_match(ismember(cur_tbt(3,:),types_mat(2,:)))==1)/sum(ismember(cur_tbt(3,:),types_mat(2,:)));
        end
    end
    full_m_tbt{m} = cur_m_tbt;
    corr_results_combine_days{m} = cur_m_corrs;
    ends_match_combine_days{m} = cur_m_ends;
end

%% Visualisation

% Combine left and right
figure
subplot(1,2,1)
types_mat = [7,10;8,11];
full_res = cell(1,2);
full_groups = cell(2,1);
groups_2 = [];
res_2 = [];
% all_violin = cell(1,2);

for m = 1:num_mice
    cur_m = corr_results_combine_days{m};
    cur_tbt = full_m_tbt{m};
    for j = 1:2
        % blc, bli, brc, bri
        full_res{j} = [full_res{j};cur_m(ismember(cur_tbt,types_mat(j,:)))];
        res_2 = [res_2;cur_m(ismember(cur_tbt,types_mat(j,:)))];
        full_groups{j} = [full_groups{j};m*ones(sum(ismember(cur_tbt,types_mat(j,:))),1)];
        groups_2 = [groups_2;j*ones(sum(ismember(cur_tbt,types_mat(j,:))),1)];
    end
end
titles_2 = ["BMI Correct","BMI Incorrect"];
    
hold on
violin_edited(full_res,'Support',[-1,1],'facecolor',[0.5,0.5,0.5],'medc',[],'edgecolor',[]);
% boxchart(groups_2,res_2,'BoxFaceColor','k','MarkerColor','k');
% title(titles_2(j));
% xlabel("Mouse")
ylim([-1,1])
xticks([1,2])
xlim([0,3])
xticklabels(["BMI Correct","BMI Incorrect"])
xtickangle(30)
ylabel("Correlations")
yticks([-1,-0.5,0,0.5,1])
yline(0,'--','LineWidth',2);
xlim([0.5,2.5])
box off

% Compare matching sign of decoder and integrated yaw at end.

% Mice Together
% figure
subplot(1,2,2)
types_vec = [7,10;8,11];
bar_mat = zeros(num_mice,2);
for m = 1:num_mice
    cur_m = ends_match_combine_days{m};
    cur_tbt = full_m_tbt{m};
    for j = 1:2
        % blc, bli, brc, bri
        bar_mat(m,j) = sum(cur_m(ismember(cur_tbt,types_vec(j,:)))==1)/sum(ismember(cur_tbt,types_vec(j,:)));
    end
end

violin_match = cell(1,2);
for m = 1:num_mice
    for d = 1:num_days
        if ~isempty(virmen_cell{m,d})
            cur_md = ends_match_cell{m,d};
            cur_tbt = tbt_cell{m,d};
            for j = 1:2
                % blc, bli, brc, bri
                violin_match{j} = [violin_match{j};sum(cur_md(ismember(cur_tbt(3,:),types_vec(j,:)))==1)/sum(ismember(cur_tbt(3,:),types_vec(j,:)))];
            end
        end
    end
end

% titles = ["BMI Correct";"BMI Incorrect"];

% hold on
% bar(mean(bar_mat),'w','LineWidth',2)
% for i = 1:2
%     plot(i*ones(num_mice,1),bar_mat(:,i),'o','Color','k','MarkerFaceColor','k','MarkerSize',marker_size)
% end

violin_edited(violin_match,'Support',[0,1+eps],'facecolor',[0.5,0.5,0.5],'medc',[],'edgecolor',[]);

% title(titles(j))
% xlabel("Mouse")
ylim([0,1])
xticks([1,2])
xlim([0,3])
xlim([0.5,2.5])
xticklabels(["BMI Correct";"BMI Incorrect"])
xtickangle(30)
ylabel("Fraction Matching")
yticks([0,0.5,1])
yline(0.5,'--','LineWidth',2);
box off

%% Hierarchical Bootstrap

% now for neuron pixel comparison
% all_p_boot = nan.*ones(4,1);
all_centres = nan.*ones(2,2);
all_sems = nan.*ones(2,2);

[all_p_boot,all_centres(1,:),all_sems(1,:)] = run_H_boot_ets(squeeze(mean_corr_results(:,:,1)), squeeze(mean_corr_results(:,:,2)),true);
[all_p_boot,all_centres(2,:),all_sems(2,:)] = run_H_boot_ets(squeeze(mean_ends_match(:,:,1)), squeeze(mean_ends_match(:,:,2)),true);

% h_boots.all_p_boot = all_p_boot;
h_boots.all_centres = all_centres;
h_boots.all_sems = all_sems;
