function [h_boots] = plot_weights_changes_corrs(va_corrs,yaw_corrs)
% 26/07/2023

% Plot results of correlating weights with neural activity changes

num_mice = size(va_corrs,2);
num_days = size(va_corrs,3);

num_sess = num_mice*num_days;
plot_off = linspace(-0.4,0.4,num_sess);

figure
scatter(plot_off+ones(1,num_sess),va_corrs(1,:),'filled','k')
hold on
plot([plot_off(1)+1,plot_off(end)+1],[mean(va_corrs(1,:),'omitnan'),mean(va_corrs(1,:),'omitnan')],'k','LineWidth',2)

mean_yaw_corrs = squeeze(mean(yaw_corrs,2,'omitnan'));
scatter(plot_off+2.*ones(1,num_sess),mean_yaw_corrs(1,:),'filled','k')
hold on
plot([plot_off(1)+2,plot_off(end)+2],[mean(mean_yaw_corrs(1,:),'omitnan'),mean(mean_yaw_corrs(1,:),'omitnan')],'k','LineWidth',2)
title(["View Angle and Ball Angular Velocity"; "Weight and Change Correlations"])
ylabel("Pearson Correlation")
yline(0,'--','LineWidth',2);
xticks([1,2])
xticklabels(["View Angle";"Ball Angular Velocity"])
ylim([-0.5,0.5])
xlim([0.5,2.5])
axis('square')

figure
scatter(plot_off+ones(1,num_sess),va_corrs(2,:)','filled','k')
hold on
plot([plot_off(1)+1,plot_off(end)+1],[mean(va_corrs(2,:),'omitnan'),mean(va_corrs(2,:),'omitnan')],'k','LineWidth',2)

mean_yaw_corrs = squeeze(mean(yaw_corrs,2,'omitnan'));
scatter(plot_off+2.*ones(1,num_sess),mean_yaw_corrs(2,:)','filled','k')
hold on
plot([plot_off(1)+2,plot_off(end)+2],[mean(mean_yaw_corrs(2,:),'omitnan'),mean(mean_yaw_corrs(2,:),'omitnan')],'k','LineWidth',2)
title(["View Angle and Ball Angular Velocity"; "Scaled Weight and Change Correlations"])
ylabel("Pearson Correlation")
yline(0,'--','LineWidth',2);
xticks([1,2])
xticklabels(["View Angle";"Ball Angular Velocity"])
ylim([-0.5,0.5])
xlim([0.5,2.5])

figure
scatter(plot_off+ones(1,num_sess),va_corrs(3,:)','filled','k')
hold on
plot([plot_off(1)+1,plot_off(end)+1],[mean(va_corrs(3,:),'omitnan'),mean(va_corrs(3,:),'omitnan')],'k','LineWidth',2)

mean_yaw_corrs = squeeze(mean(yaw_corrs,2,'omitnan'));
scatter(plot_off+2.*ones(1,num_sess),mean_yaw_corrs(3,:)','filled','k')
hold on
plot([plot_off(1)+2,plot_off(end)+2],[mean(mean_yaw_corrs(3,:),'omitnan'),mean(mean_yaw_corrs(3,:),'omitnan')],'k','LineWidth',2)
title(["View Angle and Ball Angular Velocity"; "Weight Scaling"])
ylabel("Pearson Correlation")
yline(0,'--','LineWidth',2);
xticks([1,2])
xticklabels(["View Angle";"Ball Angular Velocity"])
ylim([-0.5,0.5])
xlim([0.5,2.5])

%% Hierarchical bootstrap
% just for absolute weights with absolute changes

% get data ready
boot_samps = 1000;
num_trials = 4;

stats_data = nan.*ones(num_mice,num_days,2);
for m = 1:num_mice

    cur_data = squeeze(va_corrs(1,m,:));

    stats_data(m,1:sum(~isnan(cur_data)),1) = cur_data(~isnan(cur_data));
    
    cur_data = squeeze(mean_yaw_corrs(1,m,:));

    stats_data(m,1:sum(~isnan(cur_data)),2) = cur_data(~isnan(cur_data));

end



% p_boots = NaN(2,1);
% 
% [p_boots(1), bootstats, bootstats_center, bootstats_sem] = get_bootstrap_results_equalsamples(squeeze(stats_data(:,:,1)),squeeze(stats_data(:,:,2)),boot_samps,num_trials,'mean');
% [p_boots(2), bootstats, ~, ~] = get_bootstrap_results_equalsamples(squeeze(stats_data(:,:,2)),squeeze(stats_data(:,:,1)),boot_samps,num_trials,'mean');

[all_p_boot,all_centres,all_sems] = run_H_boot_ets(squeeze(stats_data(:,:,1)), squeeze(stats_data(:,:,2)),false);

h_boots.all_p_boot = all_p_boot;
h_boots.all_centres = all_centres;
h_boots.all_sems = all_sems;
