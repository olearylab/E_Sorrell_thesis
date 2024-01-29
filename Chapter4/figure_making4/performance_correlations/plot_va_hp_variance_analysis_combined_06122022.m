function [corr_vals,p_vals] = plot_va_hp_variance_analysis_combined_06122022(ex_mds,ex_trials,virmen_cell, hp_results_all, performance_vec, plot_inds, mice_inds, hp_cut, RMSE_all)
% 06/12/2022

% Plot figure showing:
% Example trajectories for high variance and low variance decoding.
% Correlations between open-loop variance and closed-loop variance.
% Correlations between variance and performance.
fs = 30;
mice_inds = mice_inds(plot_inds);

%% Example trajectories
ex_virmen = virmen_cell{ex_mds(1,1),ex_mds(1,2)};
ex_colours = [[0 0.4470 0.7410];[0.8500 0.3250 0.0980]];
%% high pass filter the data
hp1 = highpass(ex_virmen(7,:),hp_cut,fs);
hp2 = highpass(ex_virmen(17,:),hp_cut,fs);
%% remove invalid data
ITI = ex_virmen(8,:);
cleaned_valid = clean_valid_data(ITI);

virmen_data = ex_virmen(:,cleaned_valid);
trial_num = virmen_data(12,:);

hp1 = hp1(cleaned_valid);
hp2 = hp2(cleaned_valid);

figure
set(gcf,'position',[844,222,395,471])
ax1 = subplot(3,2,1);
% ax1 = subplot(2,4,1);
hold on
% One left trial, one right trial
% t = (1:sum(trial_num==ex_trials(1,1)))./30;
% plot(t,virmen_data(7,trial_num==ex_trials(1,1)),'--','LineWidth',2,'Color',[0 0.4470 0.7410])
% plot(t,virmen_data(17,trial_num==ex_trials(1,1)),'-','LineWidth',2,'Color',[0 0.4470 0.7410])

t = (1:sum(trial_num==ex_trials(1,2)))./30;
% plot(t,virmen_data(7,trial_num==ex_trials(1,2)),'--','LineWidth',2,'Color',[0 0.4470 0.7410])

% plot(t,virmen_data(7,trial_num==ex_trials(1,2)),'-','LineWidth',2,'Color','k')
plot(t,virmen_data(17,trial_num==ex_trials(1,2)),'-','LineWidth',2,'Color',ex_colours(1,:))
xticklabels([])
ylabel("View Angle (rad)")
title("M1")
xlim([0,10])
% legend("True","Decoded")

ax3 = subplot(3,2,3);
hold on
% plot(t,highpass(virmen_data(7,trial_num==ex_trials(1,2)),hp_cut,fs),'--','LineWidth',2,'Color',[0 0.4470 0.7410])

% plot(t,highpass(virmen_data(7,trial_num==ex_trials(1,2)),hp_cut,fs),'-','LineWidth',2,'Color','k')
plot(t,highpass(virmen_data(17,trial_num==ex_trials(1,2)),hp_cut,fs),'-','LineWidth',2,'Color',ex_colours(1,:))
%plot(t,hp1(trial_num==ex_trials(1,2)),'--','LineWidth',2,'Color',[0 0.4470 0.7410])
%plot(t,hp2(trial_num==ex_trials(1,2)),'-','LineWidth',2,'Color',[0 0.4470 0.7410])
ylabel(["1Hz High-Pass"; "View Angle (rad)"])
% xlabel("Time (s)")
xticklabels([])
xlim([0,10])
% title("True = " + var(highpass(virmen_data(7,trial_num==ex_trials(1,2)),hp_cut,fs)) + ...
%     "Decoded = " + var(highpass(virmen_data(17,trial_num==ex_trials(1,2)),hp_cut,fs)))

% Also plot true hp va variance values for each example
ax5 = subplot(3,2,5);
plot(t,highpass(virmen_data(7,trial_num==ex_trials(1,2)),hp_cut,fs),'-','LineWidth',2,'Color','k')
ylabel(["1Hz High-Pass"; "View Angle (rad)"])
xlabel("Time (s)")
xlim([0,10])
box off

%% Good mouse
ex_virmen = virmen_cell{ex_mds(2,1),ex_mds(2,2)};

%% high pass filter the data
hp1 = highpass(ex_virmen(7,:),hp_cut,fs);
hp2 = highpass(ex_virmen(17,:),hp_cut,fs);
%% remove invalid data
ITI = ex_virmen(8,:);
cleaned_valid = clean_valid_data(ITI);

virmen_data = ex_virmen(:,cleaned_valid);

hp1 = hp1(cleaned_valid);
hp2 = hp2(cleaned_valid);

trial_num = virmen_data(12,:);
ax2 = subplot(3,2,2);
hold on
% One left trial, one right trial
% t = (1:sum(trial_num==ex_trials(2,1)))./30;
% plot(t,virmen_data(7,trial_num==ex_trials(2,1)),'--','LineWidth',2,'Color','b')
% plot(t,virmen_data(17,trial_num==ex_trials(2,1)),'-','LineWidth',2,'Color','b')

t = (1:sum(trial_num==ex_trials(2,2)))./30;
% plot(t,virmen_data(7,trial_num==ex_trials(2,2)),'--','LineWidth',2,'Color',[0.4940 0.1840 0.5560])

% plot(t,virmen_data(7,trial_num==ex_trials(2,2)),'-','LineWidth',2,'Color','k')
plot(t,virmen_data(17,trial_num==ex_trials(2,2)),'-','LineWidth',2,'Color',ex_colours(2,:))
% xlabel("Time (s)")
% ylabel("View Angle (rad)")
title("M4")
xticklabels([])
yticklabels([])
linkaxes([ax1,ax2])
ax2.YAxis.TickValues = ax1.YAxis.TickValues;
xlim([0,10])

ax4 = subplot(3,2,4);
hold on
% plot(t,highpass(virmen_data(7,trial_num==ex_trials(2,2)),hp_cut,fs),'--','LineWidth',2,'Color',[0.4940 0.1840 0.5560])

% plot(t,highpass(virmen_data(7,trial_num==ex_trials(2,2)),hp_cut,fs),'-','LineWidth',2,'Color','k')
plot(t,highpass(virmen_data(17,trial_num==ex_trials(2,2)),hp_cut,fs),'-','LineWidth',2,'Color',ex_colours(2,:))
% plot(t,hp1(trial_num==ex_trials(2,2)),'--','LineWidth',2,'Color',[0.4940 0.1840 0.5560])
% plot(t,hp2(trial_num==ex_trials(2,2)),'-','LineWidth',2,'Color',[0.4940 0.1840 0.5560])
%ylabel(["High-Pass"; "View Angle (rad)"])
yticklabels([])
% xlabel("Time (s)")
xticklabels([])
ax4.YAxis.TickValues = ax3.YAxis.TickValues;
% title("True = " + var(highpass(virmen_data(7,trial_num==ex_trials(1,2)),hp_cut,fs)) + ...
    % "Decoded = " + var(highpass(virmen_data(17,trial_num==ex_trials(1,2)),hp_cut,fs)))

xlim([0,10])

ax6 = subplot(3,2,6);
plot(t,highpass(virmen_data(7,trial_num==ex_trials(2,2)),hp_cut,fs),'-','LineWidth',2,'Color','k')
% ylabel(["1Hz High-Pass"; "View Angle (rad)"])
ax6.YAxis.TickValues = ax5.YAxis.TickValues;
xlabel("Time (s)")
xlim([0,10])
box off
yticks([-0.1,0,0.1])
yticklabels([])

linkaxes([ax3,ax4,ax5,ax6])
%% Correlations
corr_vals = nan.*ones(4,1);
p_vals = nan.*ones(4,1);
% Also calculate mean true control view angle hp variance.
% Should prpbably use (95%) CIs instead of max/min
mean_true = mean(hp_results_all(plot_inds,1,1));
max_true = max(hp_results_all(plot_inds,1,1));
min_true = min(hp_results_all(plot_inds,1,1));
x_true = [min_true,min_true,max_true,max_true];
y_plot = [0.5,1,1,0.5];

%% Plotting
% performance with closed loop high passed va variance
% Change to plot different mice with different colours...
% colour_vec = [[0 0.4470 0.7410];[0 0.4470 0.7410];[0 0.4470 0.7410];[0 0.4470 0.7410];...
%    [0.8500 0.3250 0.0980];[0.8500 0.3250 0.0980];[0.8500 0.3250 0.0980];[0.8500 0.3250 0.0980];...
%    [0.9290 0.6940 0.1250];[0.9290 0.6940 0.1250];[0.9290 0.6940 0.1250];[0.9290 0.6940 0.1250];...
%    [0.4940 0.1840 0.5560];[0.4940 0.1840 0.5560];[0.4940 0.1840 0.5560];...
%    [0.4660 0.6740 0.1880];[0.4660 0.6740 0.1880];[0.4660 0.6740 0.1880];[0.4660 0.6740 0.1880]];

% colour_vec = [[0 0.4470 0.7410];[0.8500 0.3250 0.0980];[0.9290 0.6940 0.1250];[0.4940 0.1840 0.5560];[0.4660 0.6740 0.1880];[0.3010 0.7450 0.9330];[0.6350 0.0780 0.1840]];
% edited so that DW129 = M5 is green
% colour_vec = [[0 0.4470 0.7410];[0.8500 0.3250 0.0980];[0.9290 0.6940 0.1250];[0.4940 0.1840 0.5560];[0.6350 0.0780 0.1840];[0.3010 0.7450 0.9330];[0.4660 0.6740 0.1880]];
colour_vec = [[0 0 0];[0 0 0];[0 0 0];[0 0 0];[0 0 0];[0 0 0];[0 0 0]];
ex_colours = [[0 0.4470 0.7410];[0.8500 0.3250 0.0980]];
colour_vec([ex_mds(1,1),ex_mds(2,1)],:) = ex_colours;
% colour_vec = 'b';

num_dp = 3;
figure
% set(gcf,'position',[1,222,918,234])
% ax = subplot(2,4,7);
% ax = subplot(1,3,3);
scatter(hp_results_all(plot_inds,2,1),performance_vec(plot_inds,2),100,colour_vec(mice_inds,:),'filled')
hold on
[corr_vals(1),p_vals(1)] = corr(hp_results_all(plot_inds,2,1),performance_vec(plot_inds,2));
% title("\rho = " + round(cur_cor(1,2),num_dp))
ylabel(["BMI Trial"; "Fraction Correct"])
xlabel(["BMI Trial 1Hz High-Pass"; "View Angle Variance"])
% p = polyfit(hp_results_all(plot_inds,2,1),100.*performance_vec(plot_inds,2),1);
%fit_slopes(m,d,i) = p(1);
% plot(hp_results_all(plot_inds,2,1),p(2) + hp_results_all(plot_inds,2,1).*p(1),'LineWidth',2,'Color','k')
mdl = fitlm(hp_results_all(plot_inds,2,1),performance_vec(plot_inds,2));
xx = linspace(min(hp_results_all(plot_inds,2,1)),max(hp_results_all(plot_inds,2,1)),100)';
[ypred,yci] = predict(mdl,xx);
plot(xx,yci(:,1),'--','LineWidth',2,'Color','k')
plot(xx,yci(:,2),'--','LineWidth',2,'Color','k')
plot(xx,ypred,'-','LineWidth',2,'Color','k')
h = fill([xx',fliplr(xx')],[yci(:,2)',fliplr(yci(:,1)')],'k','LineStyle','none');
set(h,'facealpha',.1)

% Plot line for true view angle hp var
xline(mean_true,'LineWidth',2,'Color','k');
h = fill(x_true,y_plot,'k','LineStyle','none');
set(h,'facealpha',.1)
% xline(max_true,'--','LineWidth',2,'Color','k');
% xline(min_true,'--','LineWidth',2,'Color','k');

xlim([-0.0002,xx(end)]);
xticks([0,0.001,0.002,0.003,0.004])
ax.XAxis.Exponent = 0;
axis('square')
% closed-loop high passed va variance with open-loop hp va var
% ax = subplot(2,4,3);
figure
scatter(hp_results_all(plot_inds,1,2),hp_results_all(plot_inds,2,1),100,colour_vec(mice_inds,:),'filled')
hold on
[corr_vals(2),p_vals(2)] = corr(hp_results_all(plot_inds,1,2),hp_results_all(plot_inds,2,1));
% title("\rho = " + cur_cor(1,2))
ylabel(["Closed-Loop 1Hz High-Pass"; "View Angle Variance"])
xlabel(["Open-Loop 1Hz High-Pass"; "View Angle Variance"])
% p = polyfit(hp_results_all(plot_inds,1,2),hp_results_all(plot_inds,2,1),1);
% plot(hp_results_all(plot_inds,1,2),p(2) + hp_results_all(plot_inds,1,2).*p(1),'LineWidth',2,'Color','k')

% mdl = fitlm(hp_results_all(plot_inds,1,2),hp_results_all(plot_inds,2,1));
% xx = linspace(min(hp_results_all(plot_inds,1,2)),max(hp_results_all(plot_inds,1,2)),100)';
% [ypred,yci] = predict(mdl,xx);
% plot(xx,yci(:,1),'--','LineWidth',2,'Color','k')
% plot(xx,yci(:,2),'--','LineWidth',2,'Color','k')
% plot(xx,ypred,'-','LineWidth',2,'Color','k')

% Add y = x line
% plot(hp_results_all(plot_inds,1,2),hp_results_all(plot_inds,1,2),'-','LineWidth',2,'Color','b')

plot([0,0.0035],[0,0.0035],'--','LineWidth',2,'Color',[0.5,0.5,0.5])
xlim([0,0.0035])
ylim([0,0.0035])
xticks([0,0.001,0.002,0.003])
yticks([0,0.001,0.002,0.003])
% ax.YAxis.Exponent = 0;
ax.XAxis.Exponent = 0;
axis('square')
%% Trend with decoder accuracy.
% Plot same as for hp var, but for RMSE
% Can make input RMSE_all or R_square_all

% figure
colour_vec = [[0 0 0];[0 0 0];[0 0 0];[0 0 0];[0 0 0];[0 0 0];[0 0 0]];

% subplot(2,4,8)
% subplot(1,3,2)
figure
scatter(RMSE_all(plot_inds,2),performance_vec(plot_inds,2),100,colour_vec(mice_inds,:),'filled')
hold on
[corr_vals(3),p_vals(3)] = corr(RMSE_all(plot_inds,2),performance_vec(plot_inds,2));
% title("\rho = " + round(cur_cor(1,2),num_dp))
ylabel(["BMI Trial"; "Fraction Correct"])
xlabel(["Ball Trial"; "View Angle RMSE (rad)"]) % change if R^2
% p = polyfit(hp_results_all(plot_inds,2,1),100.*performance_vec(plot_inds,2),1);
%fit_slopes(m,d,i) = p(1);
% plot(hp_results_all(plot_inds,2,1),p(2) + hp_results_all(plot_inds,2,1).*p(1),'LineWidth',2,'Color','k')
mdl = fitlm(RMSE_all(plot_inds,2),performance_vec(plot_inds,2));
xx = linspace(min(RMSE_all(plot_inds,2)),max(RMSE_all(plot_inds,2)),100)';
[ypred,yci] = predict(mdl,xx);
plot(xx,yci(:,1),'--','LineWidth',2,'Color','k')
plot(xx,yci(:,2),'--','LineWidth',2,'Color','k')
plot(xx,ypred,'-','LineWidth',2,'Color','k')
h = fill([xx',fliplr(xx')],[yci(:,2)',fliplr(yci(:,1)')],'k','LineStyle','none');
set(h,'facealpha',.1)
xlim([xx(1),xx(end)]);

%% Trend with pitch decoder accuracy
% NB, maybe need to convert pitch into forward velocity before calculating
% RMSE for plotting... Is it ok to convert after calculating RMSE?

% subplot(2,4,4)
% subplot(1,3,1)
figure
scatter(RMSE_all(plot_inds,1),performance_vec(plot_inds,2),100,colour_vec(mice_inds,:),'filled')
hold on
[corr_vals(4),p_vals(4)] = corr(RMSE_all(plot_inds,1),performance_vec(plot_inds,2));
% title("\rho = " + round(cur_cor(1,2),num_dp))
ylabel(["BMI Trial"; "Fraction Correct"])
xlabel(["Ball Trial"; "Forward Velocity RMSE (m/s)"]) % change if R^2
% p = polyfit(hp_results_all(plot_inds,2,1),100.*performance_vec(plot_inds,2),1);
%fit_slopes(m,d,i) = p(1);
% plot(hp_results_all(plot_inds,2,1),p(2) + hp_results_all(plot_inds,2,1).*p(1),'LineWidth',2,'Color','k')
mdl = fitlm(RMSE_all(plot_inds,1),performance_vec(plot_inds,2));
xx = linspace(min(RMSE_all(plot_inds,1)),max(RMSE_all(plot_inds,1)),100)';
[ypred,yci] = predict(mdl,xx);
plot(xx,yci(:,1),'--','LineWidth',2,'Color','k')
plot(xx,yci(:,2),'--','LineWidth',2,'Color','k')
plot(xx,ypred,'-','LineWidth',2,'Color','k')
h = fill([xx',fliplr(xx')],[yci(:,2)',fliplr(yci(:,1)')],'k','LineStyle','none');
set(h,'facealpha',.1)
xlim([xx(1),xx(end)]);