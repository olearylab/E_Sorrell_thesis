function [h_boots] = plot_behav_corrs(behav_corr_cell,behav_std_corr_cell,behav_full_corr_cell)
% 25/07/2023

num_mice = size(behav_corr_cell,1);
num_days = size(behav_corr_cell,2);

corrs_mat = nan.*ones(4,4,2,num_mice,num_days);
stds_mat = nan.*ones(4,4,2,num_mice,num_days);

% Combine all into simplified 2 violin plots
v_1 = [];
v_2 = [];

for m = 1:num_mice
    for d = 1:num_days
        if~isempty(behav_corr_cell{m,d})
            for i = 1:2
                corrs_mat(:,:,i,m,d) = behav_corr_cell{m,d,i};
                stds_mat(:,:,i,m,d) = behav_std_corr_cell{m,d,i};
                v_1 = [v_1;behav_full_corr_cell{m,d,1}{1+(i-1),3+(i-1)}];
                v_2 = [v_2;behav_full_corr_cell{m,d,2}{1+(i-1),3+(i-1)}];
            end
        end
    end
end

% Plot behaviour correlation results
num_sess = num_mice*num_days;
plot_off = linspace(-0.4,0.4,num_sess);

titles = ["Forward Velocity";"Angular Velocity"];

figure
for i = 1:2
    subplot(1,2,i)
    errorbar(plot_off+ones(1,num_sess),squeeze(corrs_mat(1,1,i,:)),squeeze(stds_mat(1,1,i,:)),'o',"Color",'k')
    hold on
    plot([plot_off(1)+1,plot_off(end)+1],[mean(squeeze(corrs_mat(1,1,i,:)),'omitnan'),mean(squeeze(corrs_mat(1,1,i,:)),'omitnan')],'k','LineWidth',2)
    pp = 2;
    errorbar(plot_off+pp.*ones(1,num_sess),squeeze(corrs_mat(3,3,i,:)),squeeze(stds_mat(3,3,i,:)),'o',"Color",'k')
    hold on
    plot([plot_off(1)+pp,plot_off(end)+pp],[mean(squeeze(corrs_mat(3,3,i,:)),'omitnan'),mean(squeeze(corrs_mat(3,3,i,:)),'omitnan')],'k','LineWidth',2)
    pp = 3;
    errorbar(plot_off+pp.*ones(1,num_sess),squeeze(corrs_mat(1,3,i,:)),squeeze(stds_mat(1,3,i,:)),'o',"Color",'k')
    hold on
    plot([plot_off(1)+pp,plot_off(end)+pp],[mean(squeeze(corrs_mat(1,3,i,:)),'omitnan'),mean(squeeze(corrs_mat(1,3,i,:)),'omitnan')],'k','LineWidth',2)

    pp = 4;
    errorbar(plot_off+pp.*ones(1,num_sess),squeeze(corrs_mat(2,2,i,:)),squeeze(stds_mat(2,2,i,:)),'o',"Color",'k')
    hold on
    plot([plot_off(1)+pp,plot_off(end)+pp],[mean(squeeze(corrs_mat(2,2,i,:)),'omitnan'),mean(squeeze(corrs_mat(2,2,i,:)),'omitnan')],'k','LineWidth',2)
    pp = 5;
    errorbar(plot_off+pp.*ones(1,num_sess),squeeze(corrs_mat(4,4,i,:)),squeeze(stds_mat(4,4,i,:)),'o',"Color",'k')
    hold on
    plot([plot_off(1)+pp,plot_off(end)+pp],[mean(squeeze(corrs_mat(4,4,i,:)),'omitnan'),mean(squeeze(corrs_mat(4,4,i,:)),'omitnan')],'k','LineWidth',2)
    pp = 6;
    errorbar(plot_off+pp.*ones(1,num_sess),squeeze(corrs_mat(2,4,i,:)),squeeze(stds_mat(2,4,i,:)),'o',"Color",'k')
    hold on
    plot([plot_off(1)+pp,plot_off(end)+pp],[mean(squeeze(corrs_mat(2,4,i,:)),'omitnan'),mean(squeeze(corrs_mat(2,4,i,:)),'omitnan')],'k','LineWidth',2)
    
    
    title(titles(i))
    ylabel("Pearson Correlation")
    yline(0,'--','LineWidth',2);
    % xticks([1,2])
    % xlabel("Weights Compared")
    xticklabels([])
    axis('square')
end

%% Plot just across types comparison
% Do we want some kind of relative change?

figure
for i = 1:2
    for pp = 1:2
        errorbar(plot_off+(i-1)*2+pp.*ones(1,num_sess),squeeze(corrs_mat(1+(pp-1),3+(pp-1),i,:)),squeeze(stds_mat(1+(pp-1),3+(pp-1),i,:)),'o',"Color",'k',"MarkerFaceColor",'k')
        hold on
        plot([plot_off(1)+(i-1)*2+pp,plot_off(end)+(i-1)*2+pp],[mean(squeeze(corrs_mat(1+(pp-1),3+(pp-1),i,:)),'omitnan'),mean(squeeze(corrs_mat(1+(pp-1),3+(pp-1),i,:)),'omitnan')],'k','LineWidth',2)
    end
end
title("Correlations in behaviour across trial type")
ylabel("Pearson Correlation")
yline(0,'--','LineWidth',2);
xticks([1,2,3,4])
xticklabels(["FV Left";"FV Right";"AV Left";"AV Right"])
xtickangle(30)
% xlabel("Weights Compared")
% xticklabels([])
axis('square')
   
% combine left and right
figure
for i = 1:2
    errorbar(plot_off+i.*ones(1,num_sess),mean([squeeze(corrs_mat(1,3,i,:)),squeeze(corrs_mat(2,4,i,:))],2),mean([squeeze(stds_mat(1,3,i,:)),squeeze(stds_mat(1,3,i,:))],2),'o',"Color",'k',"MarkerFaceColor",'k')
    hold on
    plot([plot_off(1)+i,plot_off(end)+i],[mean([squeeze(corrs_mat(1,3,i,:)),squeeze(corrs_mat(2,4,i,:))],'all','omitnan'),mean([squeeze(corrs_mat(1,3,i,:)),squeeze(corrs_mat(2,4,i,:))],'all','omitnan')],'k','LineWidth',2)
end
title("Correlations in behaviour across trial type")
ylabel("Pearson Correlation")
yline(0,'--','LineWidth',2);
xticks([1,2])
xticklabels(["Forward Velocity";"Angular Velocity"])
xtickangle(30)
% xlabel("Weights Compared")
% xticklabels([])
axis('square')

% Simplify to violin
figure
v_cell{1} = v_1;
v_cell{2} = v_2;
violin(v_cell,'facecolor',[0.5,0.5,0.5],'medc',[],'edgecolor',[]);
title("Correlations in behaviour across trial type")
ylabel("Pearson Correlation")
yline(0,'--','LineWidth',2);
xticks([1,2])
xticklabels(["Forward Velocity";"Angular Velocity"])
xtickangle(30)
% xlabel("Weights Compared")
% xticklabels([])
axis('square')

% Simplify to mean and std across sessions
figure
for i = 1:2
    errorbar(i,mean([squeeze(corrs_mat(1,3,i,:)),squeeze(corrs_mat(2,4,i,:))],'all','omitnan'),std([squeeze(corrs_mat(1,3,i,:)),squeeze(corrs_mat(2,4,i,:))],0,'all','omitnan'),'o',"Color",'k',"MarkerFaceColor",'k','LineWidth',2,'MarkerSize',20)
    hold on
end
title("Correlations in behaviour across trial type")
ylabel("Pearson Correlation")
yline(0,'--','LineWidth',2);
xticks([1,2])
xlim([0.5,2.5])
xticklabels(["Forward Velocity";"Angular Velocity"])
xtickangle(30)
% xlabel("Weights Compared")
% xticklabels([])
axis('square')

% Show simplified and combined comparison
% compare within ball to wihtin bmi
figure
for i = 1:2
    errorbar(2*(i-1)+1,mean([squeeze(corrs_mat(1,1,i,:)),squeeze(corrs_mat(2,2,i,:))],'all','omitnan'),std([squeeze(corrs_mat(1,1,i,:)),squeeze(corrs_mat(2,2,i,:))],0,'all','omitnan'),'o',"Color",'k',"MarkerFaceColor",'k','LineWidth',2,'MarkerSize',20)
    hold on
end
for i = 1:2
    errorbar(2*(i-1)+2,mean([squeeze(corrs_mat(3,3,i,:)),squeeze(corrs_mat(4,4,i,:))],'all','omitnan'),std([squeeze(corrs_mat(3,3,i,:)),squeeze(corrs_mat(4,4,i,:))],0,'all','omitnan'),'o',"Color",'k',"MarkerFaceColor",'k','LineWidth',2,'MarkerSize',20)
    hold on
end
title("Correlations in behaviour within trial type")
ylabel("Pearson Correlation")
yline(0,'--','LineWidth',2);
xticks([1,2,3,4])
xlim([0.5,4.5])
xticklabels(["Ball";"BMI";"Ball";"BMI";])
xtickangle(30)
% xlabel("Weights Compared")
% xticklabels([])
axis('square')

% Show simplified and combined comparison
% compare within type to across type
figure
for i = 1:2
    errorbar(2*(i-1)+1,mean([squeeze(corrs_mat(1,1,i,:)),squeeze(corrs_mat(2,2,i,:)),squeeze(corrs_mat(3,3,i,:)),squeeze(corrs_mat(4,4,i,:))],'all','omitnan'),std([squeeze(corrs_mat(1,1,i,:)),squeeze(corrs_mat(2,2,i,:)),squeeze(corrs_mat(3,3,i,:)),squeeze(corrs_mat(4,4,i,:))],0,'all','omitnan'),'o',"Color",'k',"MarkerFaceColor",'k','LineWidth',2,'MarkerSize',20)
    hold on
end
for i = 1:2
    errorbar(2*(i-1)+2,mean([squeeze(corrs_mat(1,3,i,:)),squeeze(corrs_mat(2,4,i,:))],'all','omitnan'),std([squeeze(corrs_mat(3,3,i,:)),squeeze(corrs_mat(4,4,i,:))],0,'all','omitnan'),'o',"Color",'k',"MarkerFaceColor",'k','LineWidth',2,'MarkerSize',20)
    hold on
end
title(["Correlations in behaviour";"within vs across trial type"])
ylabel("Pearson Correlation")
yline(0,'--','LineWidth',2);
xticks([1,2,3,4])
xlim([0.5,4.5])
xticklabels(["Within";"Across";"Within";"Across"])
xtickangle(30)
% xlabel("Weights Compared")
% xticklabels([])
axis('square')

% Show simplified and combined comparison
% Plot differences from within to across



figure
for i = 1:2
    within_means = mean([squeeze(corrs_mat(1,1,i,:)),squeeze(corrs_mat(2,2,i,:)),squeeze(corrs_mat(3,3,i,:)),squeeze(corrs_mat(4,4,i,:))],2,'omitnan');
    across_means = mean([squeeze(corrs_mat(1,3,i,:)),squeeze(corrs_mat(2,4,i,:))],2,'omitnan');
    diff_means = mean(across_means-within_means,'omitnan');
    diff_stds = std(across_means-within_means,0,'omitnan');
    errorbar(i,diff_means,diff_stds,'o',"Color",'k',"MarkerFaceColor",'k','LineWidth',2,'MarkerSize',20)
    hold on
end

title(["Correlations in behaviour";"decrease across trial type"])
ylabel(["\Delta Pearson Correlation"])
% yline(0,'--','LineWidth',2);
xticks([1,2])
xlim([0.5,2.5])
xticklabels(["Forward Velocity";"Angular Velocity"])
xtickangle(30)
% xlabel("Weights Compared")
% xticklabels([])
axis('square')
box off

% plot across correlation normalised by within correlation

figure
for i = 1:2
    within_means = mean([squeeze(corrs_mat(1,1,i,:)),squeeze(corrs_mat(2,2,i,:)),squeeze(corrs_mat(3,3,i,:)),squeeze(corrs_mat(4,4,i,:))],2,'omitnan');
    across_means = mean([squeeze(corrs_mat(1,3,i,:)),squeeze(corrs_mat(2,4,i,:))],2,'omitnan');
    diff_means = mean(across_means,'omitnan')/mean(within_means,'omitnan');
    diff_stds = std(across_means,0,'omitnan')/mean(within_means,'omitnan');
    errorbar(i,diff_means,diff_stds,'o',"Color",'k',"MarkerFaceColor",'k','LineWidth',2,'MarkerSize',20)
    hold on
end

title(["Correlations in behaviour";"across trial type"])
ylabel("Normalised Pearson Correlation")
% yline(0,'--','LineWidth',2);
xticks([1,2])
xlim([0.5,2.5])
xticklabels(["Forward Velocity";"Angular Velocity"])
xtickangle(30)
% xlabel("Weights Compared")
% xticklabels([])
axis('square')
box off

%% Hierarchical bootstrap
final_means = nan.*ones(num_mice,num_days,2);
for i = 1:2
    within_means = mean([squeeze(corrs_mat(1,1,i,:,:)),squeeze(corrs_mat(2,2,i,:,:)),squeeze(corrs_mat(3,3,i,:,:)),squeeze(corrs_mat(4,4,i,:,:))],2,'omitnan');
    across_means = mean([squeeze(corrs_mat(1,3,i,:,:)),squeeze(corrs_mat(2,4,i,:,:))],2,'omitnan');
    
    cur_means = across_means-within_means;
    for m = 1:num_mice
        final_means(m,1:sum(~isnan(cur_means(m,:))),i) = cur_means(m,~isnan(cur_means(m,:)));
    end
end

% boot_samps = 1000;
% num_trials = 4;
% 
% p_boots = NaN(2,1);
% 
% [p_boots(1), bootstats, bootstats_center, bootstats_sem] = get_bootstrap_results_equalsamples(squeeze(final_means(:,:,1)),squeeze(final_means(:,:,2)),boot_samps,num_trials,'mean');
% [p_boots(2), bootstats, ~, ~] = get_bootstrap_results_equalsamples(squeeze(final_means(:,:,2)),squeeze(final_means(:,:,1)),boot_samps,num_trials,'mean');

[all_p_boot,all_centres,all_sems] = run_H_boot_ets(squeeze(final_means(:,:,1)), squeeze(final_means(:,:,2)),false);

h_boots.all_p_boot = all_p_boot;
h_boots.all_centres = all_centres;
h_boots.all_sems = all_sems;
