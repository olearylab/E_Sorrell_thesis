function [h_boots] = plot_many_corrective_decoding_check(all_res_cell)
% 26/06/2023

% Plot results of decoding check.
% 1 = train ball/ test ball. 2 = train ball/ test bmi. 3 = train bmi/ test
% ball. 4 = train bmi/ test bmi

num_mice = size(all_res_cell,1);
num_days = size(all_res_cell,2);

% What should I plot? Means? Medians? stds?
% RMSE? R^2?
% Currently mean RMSE
all_means = nan.*ones(3,8,num_mice,num_days);
all_means_R2 = nan.*ones(3,8,num_mice,num_days);
for m = 1:num_mice
    for d = 1:num_days
        if ~isempty(all_res_cell{m,d})
            
            cur_res = all_res_cell{m,d};
            
            all_means(:,:,m,d) = cur_res(:,1,:);
            
            all_means_R2(:,:,m,d) = cur_res(:,2,:);
           
        end
    end
end

%% Plot

% Slight offset of points
% Assumed 20 sessions hard coded
num_sess = num_mice*num_days;
plot_off = linspace(-0.4,0.4,num_sess);

% Yaw
% What about conversion to rad/s
circum = 64;
V = 0.32;
alpha = -50/75*circum/V;
beta = 0.05*circum/V/2.5;

figure
for i = 1:2
    scatter(plot_off+i.*ones(1,num_sess),beta.*squeeze(all_means(i,4,:))','filled','k')
    hold on
    plot([i+plot_off(1),i+plot_off(end)],beta.*[mean(squeeze(all_means(i,4,:)),'omitnan'),mean(squeeze(all_means(i,4,:)),'omitnan')],'k','LineWidth',2)
end
title(["Ball angular velocity";"decoding accuracy"])
ylabel("RMSE (rad/s)")
xlabel("Heading deviation magnitude")
xticks([1,2,3,4])
% xtickangle(30)
xticklabels(["low";"high"])
axis('square')

% Yaw, R^2, just Ball/Ball and BMI/BMI
figure
for i = 1:2
    scatter(plot_off+i.*ones(1,num_sess),squeeze(all_means_R2(i,4,:))','filled','k')
    hold on
    plot([i+plot_off(1),i+plot_off(end)],[mean(squeeze(all_means_R2(i,4,:)),'omitnan'),mean(squeeze(all_means_R2(i,4,:)),'omitnan')],'k','LineWidth',2)
end
title(["Ball angular velocity";"decoding accuracy"])
ylabel("R^2")
xticks([1,2,3,4])
% xtickangle(30)
xticklabels(["low";"high"])
xlabel("Heading deviation magnitude")
ylim([0,1])
axis('square')

%% Stats tests

%% Hierarchical bootstrap
means_ready = nan.*ones(2,num_mice,num_days);
for m = 1:num_mice
    cur_means = squeeze(all_means_R2(1:2,4,m,:));
    means_ready(:,m,1:sum(~isnan(cur_means(1,:)))) = cur_means(:,~isnan(cur_means(1,:)));
end

% boot_samps = 1000;
% num_trials = 4;
% 
% p_boots = NaN(2,1);
% 
% [p_boots(1), bootstats, bootstats_center, bootstats_sem] = get_bootstrap_results_equalsamples(squeeze(means_ready(1,:,:)),squeeze(means_ready(2,:,:)),boot_samps,num_trials,'mean');
% [p_boots(2), bootstats, ~, ~] = get_bootstrap_results_equalsamples(squeeze(means_ready(2,:,:)),squeeze(means_ready(1,:,:)),boot_samps,num_trials,'mean');

[all_p_boot,all_centres,all_sems] = run_H_boot_ets(squeeze(means_ready(1,:,:)), squeeze(means_ready(2,:,:)),false);

h_boots.all_p_boot = all_p_boot;
h_boots.all_centres = all_centres;
h_boots.all_sems = all_sems;
