function [hp_results_all,yaw_results_all,performance_vec] = plot_variability_performance_correlations(virmen_cell,tbt_cell,summary_cell,mice_vec,yaw_offsets,hp_filter_cut,var_samps)
% 16/05/2022

% Function for checking correlations between yaw variances/view angle
% fluctuations etc. with behavioural performance.

% cells are 1d
num_sessions = length(virmen_cell);

hp_results_all = nan.*ones(num_sessions,2,2);
yaw_results_all = nan.*ones(num_sessions,2,2);
% Keep total average, maybe split?
performance_vec = nan.*ones(num_sessions,2);

%%
for n = 1:num_sessions
    virmen_data = virmen_cell{n};
    tbt_details = tbt_cell{n};
    tbt_summary = summary_cell{n};
    
    % [all_processed,all_vars,final_results] = trial_averaged_variance_measures(virmen_cell{m},tbt_cell{m},31,hp_filter_cut,31);
    [all_processed,all_vars,final_results,final_yaw_results] = trial_averaged_variance_measures_full(virmen_data,tbt_details,var_samps,hp_filter_cut,var_samps,yaw_offsets(n));
    hp_results_all(n,:,:) = final_results(2,:,:);
    yaw_results_all(n,:,:) = final_yaw_results; % control/bmi x observed/decoded
    performance_vec(n,1) = sum(tbt_summary([1,4]))./sum(tbt_summary(1:6));
    performance_vec(n,2) = sum(tbt_summary([7,10]))./sum(tbt_summary(7:end));
end

%% Plotting
% performance with closed loop high passed va variance
figure
scatter(hp_results_all(:,2,1),100.*performance_vec(:,2),'filled')
hold on
cur_cor = corrcoef(hp_results_all(:,2,1),performance_vec(:,2));
title("Correlation Value = " + cur_cor)
ylabel("BMI Performance")
xlabel("Closed-Loop High-Pass View Angle Variance")
p = polyfit(hp_results_all(:,2,1),100.*performance_vec(:,2),1);
%fit_slopes(m,d,i) = p(1);
plot(hp_results_all(:,2,1),p(2) + hp_results_all(:,2,1).*p(1),'LineWidth',2,'Color','k')

% closed-loop high passed va variance with open-loop hp va var
figure
scatter(hp_results_all(:,1,2),hp_results_all(:,2,1),'filled')
hold on
cur_cor = corrcoef(hp_results_all(:,1,2),hp_results_all(:,2,1));
title("Correlation Value = " + cur_cor)
ylabel("Open-Loop High-Pass View Angle Variance")
xlabel("Closed-Loop High-Pass View Angle Variance")
p = polyfit(hp_results_all(:,1,2),hp_results_all(:,2,1),1);
plot(hp_results_all(:,1,2),p(2) + hp_results_all(:,1,2).*p(1),'LineWidth',2,'Color','k')

% performance with open-loop va var
figure
scatter(hp_results_all(:,1,2),100.*performance_vec(:,2),'filled')
hold on
cur_cor = corrcoef(hp_results_all(:,1,2),performance_vec(:,2));
title("Correlation Value = " + cur_cor)
ylabel("BMI Performance")
xlabel("Open-Loop High-Pass View Angle Variance")
p = polyfit(hp_results_all(:,1,2),100.*performance_vec(:,2),1);
plot(hp_results_all(:,1,2),p(2) + hp_results_all(:,1,2).*p(1),'LineWidth',2,'Color','k')

% performance with closed-loop yaw movingvar
figure
scatter(yaw_results_all(:,2,2),100.*performance_vec(:,2),'filled')
hold on
cur_cor = corrcoef(yaw_results_all(:,2,2),performance_vec(:,2));
title("Correlation Value = " + cur_cor)
ylabel("BMI Performance")
xlabel("Closed-Loop Decoded Yaw Moving Variance")
p = polyfit(yaw_results_all(:,2,2),100.*performance_vec(:,2),1);
plot(yaw_results_all(:,2,2),p(2) + yaw_results_all(:,2,2).*p(1),'LineWidth',2,'Color','k')

% open-loop with closed-loop yaw movingvar
figure
scatter(yaw_results_all(:,1,2),yaw_results_all(:,2,2),'filled')
hold on
cur_cor = corrcoef(yaw_results_all(:,1,2),yaw_results_all(:,2,2));
title("Correlation Value = " + cur_cor)
ylabel("Closed-Loop Decoded Yaw Moving Variance")
xlabel("Open-Loop Decoded Yaw Moving Variance")
p = polyfit(yaw_results_all(:,1,2),yaw_results_all(:,2,2),1);
plot(yaw_results_all(:,1,2),p(2) + yaw_results_all(:,1,2).*p(1),'LineWidth',2,'Color','k')