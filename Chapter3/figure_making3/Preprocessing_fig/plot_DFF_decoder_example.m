function [] = plot_DFF_decoder_example(results_struct_raw,results_struct_dff)
% 25/05/2023

% 2x2 figure of true and decoder output at start, and near end of testing
% session

% colours_vec = [[0 0.4470 0.7410];[0.8500 0.3250 0.0980];[0.9290 0.6940 0.1250];[0.4940 0.1840 0.5560];[0.4660 0.6740 0.1880]];
line_width = 2;
t = (1:size(results_struct_raw.xtest(results_struct_raw.test_valid,:),1))/results_struct_raw.model_params.fs;

ax1 = subplot(2,2,1);

max1 = max(results_struct_raw.xtest(:,1));
max2 = max(results_struct_raw.xprediction(:,1));

min1 = min(results_struct_raw.xtest(:,1));
min2 = min(results_struct_raw.xprediction(:,1));

baseline = min([min1,min2]);
max_val = max([max1,max2]);

hold on
plot(t,0.01.*results_struct_raw.xtest(results_struct_raw.test_valid,1),'LineWidth',line_width,'Color',[0.5,0.5,0.5])

plot(t,0.01.*results_struct_raw.xprediction(results_struct_raw.test_valid,1),'LineWidth',line_width,'Color',[0.4940 0.1840 0.5560])
title("Start of Testing")
ylabel("Y Position (m)")
% xlabel("Time (s)")
xlim([0,60])
% ylim([baseline,max_val]);

xticks([0, 60])
xticklabels({'0','60'})
% yticks(ytick_vec(i,:))
%yticklabels({'-1','0','1'})
% ylabel(ylabel_vec(i))
ylim([-1.20,5.10])

box off

ax2 = subplot(2,2,2);

hold on
plot(t,0.01.*results_struct_raw.xtest(results_struct_raw.test_valid,1),'LineWidth',line_width,'Color',[0.5,0.5,0.5])

plot(t,0.01.*results_struct_raw.xprediction(results_struct_raw.test_valid,1),'LineWidth',line_width,'Color',[0.4940 0.1840 0.5560])
title("End of Testing")
% ylabel("Y Position (m)")
% xlabel("Time (s)")
xlim([t(end)-60, t(end)])
% ylim([baseline,max_val]);

% xticks([0, 60])
% xticklabels({'0','60'})
% yticks(ytick_vec(i,:))
%yticklabels({'-1','0','1'})
% ylabel(ylabel_vec(i))
ylim([-1.20,5.10])

ax3 = subplot(2,2,3);

hold on
plot(t,0.01.*results_struct_dff.xtest(results_struct_dff.test_valid,1),'LineWidth',line_width,'Color',[0.5,0.5,0.5])

plot(t,0.01.*results_struct_dff.xprediction(results_struct_dff.test_valid,1),'LineWidth',line_width,'Color',[0.4940 0.1840 0.5560])
title("Start of Testing")
ylabel("Y Position (m)")
xlabel("Time (s)")
xlim([0,60])
% ylim([baseline,max_val]);

xticks([0, 60])
xticklabels({'0','60'})
% yticks(ytick_vec(i,:))
%yticklabels({'-1','0','1'})
% ylabel(ylabel_vec(i))
ylim([-1.20,5.10])

ax4 = subplot(2,2,4);

hold on
plot(t,0.01.*results_struct_dff.xtest(results_struct_dff.test_valid,1),'LineWidth',line_width,'Color',[0.5,0.5,0.5])

plot(t,0.01.*results_struct_dff.xprediction(results_struct_dff.test_valid,1),'LineWidth',line_width,'Color',[0.4940 0.1840 0.5560])
title("End of Testing")
% ylabel("Y Position (m)")
xlabel("Time (s)")
xlim([t(end)-60, t(end)])
% ylim([baseline,max_val]);
ylim([-1.20,5.10])

% xticks([0, 60])
% xticklabels({'0','60'})
% yticks(ytick_vec(i,:))
%yticklabels({'-1','0','1'})
% ylabel(ylabel_vec(i))

linkaxes([ax1,ax3],'x')
linkaxes([ax2,ax4],'x')

