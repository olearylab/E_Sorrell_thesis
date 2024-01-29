function [h_boots] = position_decoding_summary_plot_09012023(results_struct_cell, virmen_cell, tbt_cell, ind)

% Plotting of position decoding
% TODO: Allow for exclusion of specific days.

num_mice = size(results_struct_cell,1);
num_days = size(results_struct_cell,2);
pos_scale = 0.74; % Not sure this works for standard deviation. Need to check my 3F3.

%% Offline Position Decoding
% get trial by trial accuracy
all_trials_cell = cell(num_mice,num_days);
for m = 1:num_mice
    for d = 1:num_days
        if ~isempty(virmen_cell{m,d})
            [all_trials_res,trials_summary_means,trials_summary_stds] = trial_by_trial_accuracy(results_struct_cell{m,d}, virmen_cell{m,d}, tbt_cell{m,d});
            all_trials_cell{m,d} = all_trials_res;
        end
    end
end

% Need trials sorted according to trial type we want to plot by. And
% combine across days, and across all mice.
% Will just work with RMSE for now, potentially R^2 as well.
full_m_trials = cell(num_mice,1);
full_all_trials = [];
full_m_tbt = cell(num_mice,1);
full_all_tbt = [];
trials_means = NaN(num_mice,num_days,6);
types_mat = [1,4;7,10;2,5;8,11;3,6;9,12];
for m = 1:num_mice
    cur_m_trials = [];
    cur_m_tbt = [];
    for d = 1:num_days
        if ~isempty(all_trials_cell{m,d})
            cur_tbt = tbt_cell{m,d};
            cur_trials = all_trials_cell{m,d};
            cur_m_trials = [cur_m_trials;cur_trials];
            full_all_trials = [full_all_trials;cur_trials];

            for i = 1:6
                trials_means(m,d,i) = mean(cur_trials(ismember(cur_tbt(3,:),types_mat(i,:)),1,ind),'omitnan')*pos_scale;
            end
            % Create new tbt_vector. For now not distinguishing direction
            % (potentially come back and also do this for bias checking).
            for i = 1:3
                cur_tbt(3,(ismember(cur_tbt(3,:),[1,4]+(i-1)))) = i;
            end
            for i = 1:3
                cur_tbt(3,(ismember(cur_tbt(3,:),[7,10]+(i-1)))) = i+3;
            end

            cur_m_tbt = [cur_m_tbt,cur_tbt(3,:)];
            full_all_tbt = [full_all_tbt,cur_tbt(3,:)];
        end
    end
    full_m_trials{m} = cur_m_trials;
    full_m_tbt{m} = cur_m_tbt;
end

% Plot box and whiskers for each mouse separately, showing normal correct,
% normal incorrect, normal timeout, BMIc, BMIi, BMIt. Summarise across all
% days. Potentially plot faded mean and standard deviation for each day as well. 

% figure
% for m = 1:num_mice
%     cur_res = full_m_trials{m};
%     subplot(1,num_mice,m)
%     hold on
%     boxplot(squeeze(cur_res(:,1,ind))*pos_scale,full_m_tbt{m})
%     title("Mouse " + m)
%     if m ==1
%         ylabel("RMSE (cm)")
%     else
%         yticks([])
%     end
%     ylim([0,200*pos_scale])
%     xlabel("Trial Type")
%     box off
%     
% end
% 
% % Now with individual day additions
% figure
% for m = 1:num_mice
%     cur_res = full_m_trials{m};
%     subplot(1,num_mice,m)
%     hold on
%     boxplot(squeeze(cur_res(:,1,ind))*pos_scale,full_m_tbt{m})
%     title("Mouse " + m)
%     if m ==1
%         ylabel("RMSE (cm)")
%     else
%         yticks([])
%     end
%     ylim([0,200*pos_scale])
%     xlabel("Trial Type")
%     box off
%     
%     for d = 1:num_days
%         cur_res = all_trials_cell{m,d};
%         cur_tbt = tbt_cell{m,d};
%         for i = 1:3
%             cur_tbt(3,(ismember(cur_tbt(3,:),[1,4]+(i-1)))) = i;
%         end
%         for i = 1:3
%             cur_tbt(3,(ismember(cur_tbt(3,:),[7,10]+(i-1)))) = i+3;
%         end
%         
%         for i = 1:6
%             errorbar(i,pos_scale*mean(squeeze(cur_res(cur_tbt(3,:)==i,1,ind))),std(pos_scale*squeeze(cur_res(cur_tbt(3,:)==i,1,ind))),'o','Color','k','MarkerFaceColor','k')
%         end
%         
%     end
%     
% end

% Plot dot and error bars with boxplot, each dot is single mouse results.
% figure
% hold on
% boxplot(squeeze(full_all_trials(:,1,ind))*pos_scale,full_all_tbt)
% title("Position Decoding Accuracy")
% ylabel("RMSE (cm)")
% ylim([0,200*pos_scale])
% xlabel("Trial Type")
% box off
% all_means = zeros(num_mice,6);
% for m = 1:num_mice
%     cur_res = full_m_trials{m};
%     cur_tbt = full_m_tbt{m};
%     for i = 1:6
%         errorbar(i,pos_scale*mean(squeeze(cur_res(cur_tbt==i,1,ind))),std(pos_scale*squeeze(cur_res(cur_tbt==i,1,ind))),'o','Color','k','MarkerFaceColor','k')
%         all_means(m,i) = pos_scale*mean(squeeze(cur_res(cur_tbt==i,1,ind)));
%     end  
% end
% xticks([1,2,3,4,5,6])
% xticklabels(["NC";"NI";"NT";"BC";"BI";"BT"])
% 
% % Alternatively, plot bar of means instead of overall boxplots. 
% figure
% hold on
% % boxplot(squeeze(full_all_trials(:,1,ind))*pos_scale,full_all_tbt)
% title("Position Decoding Accuracy")
% ylabel("RMSE (cm)")
% ylim([0,200*pos_scale])
% xlabel("Trial Type")
% box off
% bar(mean(all_means,'omitnan'),'w','LineWidth',2)
% for m = 1:num_mice
%     cur_res = full_m_trials{m};
%     cur_tbt = full_m_tbt{m};
%     for i = 1:6
%         errorbar(i,pos_scale*mean(squeeze(cur_res(cur_tbt==i,1,ind))),std(pos_scale*squeeze(cur_res(cur_tbt==i,1,ind))),'o','Color','k','MarkerFaceColor','k')
%         % all_means(m,i) = pos_scale*mean(squeeze(cur_res(cur_tbt==i,1,ind)));
%     end  
% end
% xticks([1,2,3,4,5,6])
% xticklabels(["NC";"NI";"NT";"BC";"BI";"BT"])

% Add statistical significance checking.

%% Violin plot version

violin_cell = cell(1,6);
for i = 1:6
    violin_cell{i} = squeeze(full_all_trials(full_all_tbt==i,1,ind))*pos_scale;
end

off_size = 0.1;

m_offset = -off_size*(num_mice-1)/2:off_size:off_size*(num_mice-1)/2;

box_pos = [1,2];

figure
hold on
violin(violin_cell,'facecolor',[0.5,0.5,0.5],'medc',[],'edgecolor',[]);
title("Position Decoding Accuracy")
ylabel("RMSE (cm)")
ylim([0,200*pos_scale])
xlabel("Trial Type")
box off
all_means = zeros(num_mice,6);
for m = 1:num_mice
    cur_res = full_m_trials{m};
    cur_tbt = full_m_tbt{m};
    for i = 1:6
        % errorbar(i+m_offset(m),pos_scale*mean(squeeze(cur_res(cur_tbt==i,1,ind))),std(pos_scale*squeeze(cur_res(cur_tbt==i,1,ind))),'o','Color','k','MarkerFaceColor','k','LineWidth',1,'MarkerSize',10)
        all_means(m,i) = pos_scale*mean(squeeze(cur_res(cur_tbt==i,1,ind)));
    end  
end
xticks([1,2,3,4,5,6])
xticklabels(["NC";"NI";"NT";"BC";"BI";"BT"])

v_alt = cell(1,6);
v_alt{1} = violin_cell{1};
v_alt{2} = violin_cell{4};
v_alt{3} = violin_cell{2};
v_alt{4} = violin_cell{5};
v_alt{5} = violin_cell{3};
v_alt{6} = violin_cell{6};

figure
hold on
violin(v_alt,'facecolor',[0.5,0.5,0.5],'medc',[],'edgecolor',[]);
title("Position Decoding Accuracy")
ylabel("RMSE (cm)")
ylim([0,200*pos_scale])
xlabel("Trial Type")
box off
all_means = zeros(num_mice,6);
for m = 1:num_mice
    cur_res = full_m_trials{m};
    cur_tbt = full_m_tbt{m};
    for i = 1:6
        % errorbar(i+m_offset(m),pos_scale*mean(squeeze(cur_res(cur_tbt==i,1,ind))),std(pos_scale*squeeze(cur_res(cur_tbt==i,1,ind))),'o','Color','k','MarkerFaceColor','k','LineWidth',1,'MarkerSize',10)
        all_means(m,i) = pos_scale*mean(squeeze(cur_res(cur_tbt==i,1,ind)));
    end  
end
xticks([1,2,3,4,5,6])
xticklabels(["Ball";"BMI";"Ball";"BMI";"Ball";"BMI"])
axis('square')
%% Hierarchical Bootstrap

all_p_boot = nan.*ones(2,1);
all_centres = nan.*ones(2,2);
all_sems = nan.*ones(2,2);

for i = 1:2
    [all_p_boot(i),all_centres(i,:),all_sems(i,:)] = run_H_boot_ets(squeeze(trials_means(:,:,2*(i-1)+1)), squeeze(trials_means(:,:,2*(i-1)+2)),true);
end

h_boots.all_p_boot = all_p_boot;
h_boots.all_centres = all_centres;
h_boots.all_sems = all_sems;
