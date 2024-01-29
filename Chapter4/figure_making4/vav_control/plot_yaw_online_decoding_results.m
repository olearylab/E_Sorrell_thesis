function [sig_results_tot,boot_stats] = plot_yaw_online_decoding_results(t_types_summary_train,t_types_summary_cell)
% 12/04/2023

% Plot single mouse results of yaw decoding BMI

num_days = length(t_types_summary_cell);

final_res = nan.*ones(4,5);
for i = 1:2
    final_res(i,1) = t_types_summary_train(3*(i-1)+1)/sum(t_types_summary_train(3*(i-1)+1:3*i));
end

for d = 1:num_days
    cur_res = t_types_summary_cell{d};
    for i = 1:4
        final_res(i,d+1) = cur_res(3*(i-1)+1)/sum(cur_res(3*(i-1)+1:3*i));
    end
end

figure
ball_colour = [0.4940 0.1840 0.5560];
bmi_colour = [0.4660 0.6740 0.1880];
% Will need to edit as I go to do exactly what I want.

set(0,'DefaultAxesFontSize',20)
lines_vec = ["-o";"--^";"-o";"--^"];
for i = 1:4
    if i < 3
            plot(1:4,final_res(i,2:end)',lines_vec(i),'LineWidth',2,'Color',ball_colour,'MarkerFaceColor',ball_colour,'MarkerSize',10)
            hold on
    else
            plot(1:4,final_res(i,2:end)',lines_vec(i),'LineWidth',2,'Color',bmi_colour,'MarkerFaceColor',bmi_colour,'MarkerSize',10)  
    end
end
for i = 1:2
    plot(0,final_res(i,1)',lines_vec(i),'LineWidth',2,'Color',ball_colour,'MarkerFaceColor',ball_colour,'MarkerSize',10)
end
yline(0.5,'LineWidth',2);
title("View Angle Velocity BMI Results")
yticks([0, 0.25, 0.5, 0.75, 1])
ylabel("Fraction Correct")
% yticks([0, 0.25, 0.5, 0.75, 1])
yticklabels([0, 0.25, 0.5, 0.75, 1])
% Some hardcoded stuff for labels etc.
xlabel("Day")
xlim([-0.5,4+0.5])
ylim([0,1.1])
xticks(0:4)
xticklabels({'Train','1','2','3','4'})
box off
legend('Ball Left', 'Ball Right', 'BMI Left', 'BMI Right')

%%

summary_results = nan.*ones(1,length(t_types_summary_cell),2);
tot_b = nan.*ones(length(t_types_summary_cell),2);
tot_n = nan.*ones(length(t_types_summary_cell),2);
for d = 1:length(t_types_summary_cell)
    cur_res = t_types_summary_cell{d};
    summary_results(1,d,1) = sum(cur_res([1,4]),'all')/sum(cur_res(1:6),'all');
    summary_results(1,d,2) = sum(cur_res([7,10]),'all')/sum(cur_res(7:end),'all');

    tot_b(d,1) = sum(cur_res([7,10]),'all');
    tot_b(d,2) = sum(cur_res(7:end),'all');
    tot_n(d,1) = sum(cur_res([1,4]),'all');
    tot_n(d,2) = sum(cur_res(1:6),'all');
end

% sig_results = nan.*ones(length(trial_types_summary_cell),num_days-1,2);
% sig_results_train = nan.*ones(length(trial_types_summary_cell),1);
% summed sig check
tot_b_all = squeeze(sum(tot_b,1,'omitnan'));
tot_n_all = squeeze(sum(tot_n,1,'omitnan'));

sig_results_tot = nan.*ones(2,1);

sig_results_tot(1) = myBinomTest(tot_n_all(1),tot_n_all(2),0.5,'one');
sig_results_tot(2) = myBinomTest(tot_b_all(1),tot_b_all(2),0.5,'one');


%% Hierarchical bootstrap
% nans already is correct place
% traing only has one day? regular bootstrap?
rng(1)
boot_samps = 1000;
num_trials = 4;

bootstats = get_bootstrapped_equalsamples(squeeze(summary_results(:,:,1)),boot_samps,num_trials,'mean');

boot_stats.ball_center = mean(bootstats);
boot_stats.ball_std = std(bootstats);

bootstats = get_bootstrapped_equalsamples(squeeze(summary_results(:,:,2)),boot_samps,num_trials,'mean');

boot_stats.bmi_center = mean(bootstats);
boot_stats.bmi_std = std(bootstats);