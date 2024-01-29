function [stats_tests,h_boots] = plot_many_decoding_check(all_res_cell)
% 26/06/2023

% Plot results of decoding check.
% 1 = train ball/ test ball. 2 = train ball/ test bmi. 3 = train bmi/ test
% ball. 4 = train bmi/ test bmi

num_mice = size(all_res_cell,1);
num_days = size(all_res_cell,2);

% What should I plot? Means? Medians? stds?
% RMSE? R^2?
% Currently mean RMSE
all_means = nan.*ones(4,8,num_mice,num_days);
all_means_R2 = nan.*ones(4,8,num_mice,num_days);
for m = 1:num_mice
    for d = 1:num_days
        if ~isempty(all_res_cell{m,d})
            
            cur_res = all_res_cell{m,d};
            
            all_means(:,:,m,d) = mean(cur_res(:,:,1,:),2,'omitnan');
            
            all_means_R2(:,:,m,d) = mean(cur_res(:,:,2,:),2,'omitnan');
           
        end
    end
end

%% Plot

% Slight offset of points
% Assumed 20 sessions hard coded
num_sess = num_mice*num_days;
plot_off = linspace(-0.4,0.4,num_sess);

% Reorder: 1. ball/ball. 2. bmi/bmi. 3. ball/bmi. 4. bmi/ball
all_means_orig = all_means;
all_means(2,:,:,:) = all_means_orig(4,:,:,:);
all_means(3,:,:,:) = all_means_orig(2,:,:,:);
all_means(4,:,:,:) = all_means_orig(3,:,:,:);

all_means_orig_R2 = all_means_R2;
all_means_R2(2,:,:,:) = all_means_orig_R2(4,:,:,:);
all_means_R2(3,:,:,:) = all_means_orig_R2(2,:,:,:);
all_means_R2(4,:,:,:) = all_means_orig_R2(3,:,:,:);

% Yposition
figure
for i = 1:4
    scatter(plot_off+i.*ones(1,num_sess),squeeze(all_means(i,1,:))','filled','k')
    hold on
    plot([i+plot_off(1),i+plot_off(end)],[mean(squeeze(all_means(i,1,:)),'omitnan'),mean(squeeze(all_means(i,1,:)),'omitnan')],'k','LineWidth',2)
end
title("Y Position")
ylabel("RMSE (vu)")
xticks([1,2,3,4])
xtickangle(30)
xticklabels(["Ball/Ball";"BMI/BMI";"Ball/BMI";"BMI/Ball"])

% Pitch
figure
for i = 1:4
    scatter(plot_off+i.*ones(1,num_sess),squeeze(all_means(i,3,:))','filled','k')
    hold on
    plot([i+plot_off(1),i+plot_off(end)],[mean(squeeze(all_means(i,3,:)),'omitnan'),mean(squeeze(all_means(i,3,:)),'omitnan')],'k','LineWidth',2)
end
title("Pitch")
ylabel("RMSE (V)")
xticks([1,2,3,4])
xtickangle(30)
xticklabels(["Ball/Ball";"BMI/BMI";"Ball/BMI";"BMI/Ball"])

% Yaw
% What about conversion to rad/s
circum = 64;
V = 0.32;
alpha = -50/75*circum/V;
beta = 0.05*circum/V/2.5;

figure
for i = 1:4
    scatter(plot_off+i.*ones(1,num_sess),beta.*squeeze(all_means(i,4,:))','filled','k')
    hold on
    plot([i+plot_off(1),i+plot_off(end)],beta.*[mean(squeeze(all_means(i,4,:)),'omitnan'),mean(squeeze(all_means(i,4,:)),'omitnan')],'k','LineWidth',2)
end
title(["Ball Angular Velocity";"Decoding Accuracy"])
ylabel("RMSE (rad/s)")
xticks([1,2,3,4])
xtickangle(30)
xticklabels(["Ball/Ball";"BMI/BMI";"Ball/BMI";"BMI/Ball"])
axis('square')

% Yaw, R^2, just Ball/Ball and BMI/BMI
figure
for i = 1:2
    scatter(plot_off+i.*ones(1,num_sess),squeeze(all_means_R2(i,4,:))','filled','k')
    hold on
    plot([i+plot_off(1),i+plot_off(end)],[mean(squeeze(all_means_R2(i,4,:)),'omitnan'),mean(squeeze(all_means_R2(i,4,:)),'omitnan')],'k','LineWidth',2)
end
title(["Ball Angular Velocity";"Decoding Accuracy"])
ylabel("R^2")
xticks([1,2,3,4])
% xtickangle(30)
xticklabels(["Ball";"BMI"])
xlabel("Trial Type")
ylim([0,1])
axis('square')

%% Stats tests
%% Statistical tests
% Perform two sided-Wilcoxon signed rank test, and two sided paired sample
% t-test, using each session as a data point (19 sessions). Also for each
% mouse as data point
all_means = beta.*all_means;
p_wilc = zeros(4,4);
% pr_wilc = zeros(2,4);

p_t = zeros(4,4);
% pr_t = zeros(2,4);

p_wilc_all = zeros(4,4);
p_t_all = zeros(4,4);

p_means = zeros(1,4);
p_means_all = zeros(1,4);
p_medians = zeros(1,4);
p_medians_all = zeros(1,4);

all_m_means = squeeze(mean(all_means,4,'omitnan'));

for i = 1:4
    for j = 1:4
        p_wilc(i,j) = signrank(squeeze(all_m_means(i,4,:)),squeeze(all_m_means(j,4,:)));

        p_wilc_all(i,j) = signrank(squeeze(all_means(i,4,:)),squeeze(all_means(j,4,:)));
                  
        [~,p_t(i,j)] = ttest(squeeze(all_m_means(i,4,:)),squeeze(all_m_means(j,4,:)));
        [~,p_t_all(i,j)] = ttest(squeeze(all_means(i,4,:)),squeeze(all_means(j,4,:)));
    end
    p_means(i) = mean(squeeze(all_m_means(i,4,:)),'omitnan');
    p_means_all(i) = mean(squeeze(all_means(i,4,:)),'omitnan');
    p_medians(i) = median(squeeze(all_m_means(i,4,:)),'omitnan');
    p_medians_all(i) = median(squeeze(all_means(i,4,:)),'omitnan');

end
stats_tests.p_wilc = p_wilc;
stats_tests.p_wilc_all = p_wilc_all;
stats_tests.p_means_all = p_means_all;
stats_tests.p_medians_all = p_medians_all;
stats_tests.p_means = p_means;
stats_tests.p_medians = p_medians;
stats_tests.p_t = p_t;
stats_tests.p_t_all = p_t_all;

% Display stats test results and medians
% Display results for R^2 at the moment
row_names = ["Ball/Ball";"BMI/BMI";"Ball/BMI";"BMI/Ball"];
dim_names = ["Median All";"Mean All";"Median Mice";"Mean Mice"];
T = table(stats_tests.p_medians_all',stats_tests.p_means_all',stats_tests.p_medians',stats_tests.p_means','RowNames',row_names,'VariableNames',dim_names);
disp(T) 

row_names = ["Ball/Ball";"BMI/BMI";"Ball/BMI";"BMI/Ball"];
dim_names = ["Ball/Ball";"BMI/BMI";"Ball/BMI";"BMI/Ball"];
T = table(stats_tests.p_wilc_all(:,1),stats_tests.p_wilc_all(:,2),stats_tests.p_wilc_all(:,3),stats_tests.p_wilc_all(:,4),'RowNames',row_names,'VariableNames',dim_names);
disp(T)

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
