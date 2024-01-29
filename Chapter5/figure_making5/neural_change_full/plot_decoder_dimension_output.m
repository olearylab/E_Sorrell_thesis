function [h_boots] = plot_decoder_dimension_output(all_means_cell,all_stds_cell)
% 27/07/2023

% plot mean absolute activity through decoder, and through absolute decoder

num_mice = size(all_means_cell,1);
num_days = size(all_means_cell,2);

res_mat = nan.*ones(4,2,num_mice,num_days);

std_mat = nan.*ones(4,2,num_mice,num_days);

for m = 1:num_mice
    for d = 1:num_days
        if ~isempty(all_means_cell{m,d})
            res_mat(:,:,m,d) = all_means_cell{m,d};
            std_mat(:,:,m,d) = all_stds_cell{m,d};
        end
    end
end

% res_mat = res_mat(:,:,:);
res_mat_on = res_mat([1,3],:,:,:);
res_mat_off = res_mat([2,4],:,:,:);

std_mat_on = std_mat([1,3],:,:,:);
std_mat_off = std_mat([2,4],:,:,:);

num_sess = num_mice*num_days;
plot_off = linspace(-0.4,0.4,num_sess);

figure
for i = 1:2
    for j = 1:2
        subplot(1,2,1)
        scatter(plot_off+((i-1)*2+j)*ones(1,num_sess),res_mat_on(i,j,:),'filled','k')
        hold on
        plot([plot_off(1)+((i-1)*2+j),plot_off(end)+((i-1)*2+j)],[mean(res_mat_on(i,j,:),'omitnan'),mean(res_mat_on(i,j,:),'omitnan')],'k','LineWidth',2)  
        
        subplot(1,2,2)
        scatter(plot_off+((i-1)*2+j)*ones(1,num_sess),res_mat_off(i,j,:),'filled','k')
        hold on
        plot([plot_off(1)+((i-1)*2+j),plot_off(end)+((i-1)*2+j)],[mean(res_mat_off(i,j,:),'omitnan'),mean(res_mat_off(i,j,:),'omitnan')],'k','LineWidth',2) 
    end
end

figure

for i = 1:2
    subplot(1,2,1)
    scatter(plot_off+(i)*ones(1,num_sess),(res_mat_on(i,2,:)-res_mat_on(i,1,:)),'filled','k')
    hold on
    plot([plot_off(1)+(i),plot_off(end)+(i)],[mean((res_mat_on(i,2,:)-res_mat_on(i,1,:)),'omitnan'),mean((res_mat_on(i,2,:)-res_mat_on(i,1,:)),'omitnan')],'k','LineWidth',2)  

    subplot(1,2,2)
    scatter(plot_off+(i)*ones(1,num_sess),(res_mat_off(i,2,:)-res_mat_off(i,1,:)),'filled','k')
    hold on
    plot([plot_off(1)+(i),plot_off(end)+(i)],[mean((res_mat_off(i,2,:)-res_mat_off(i,1,:)),'omitnan'),mean((res_mat_off(i,2,:)-res_mat_off(i,1,:)),'omitnan')],'k','LineWidth',2) 
end

subplot(1,2,1)
yline(0,'--','LineWidth',2);
title("View Angle")
ylabel(["Mean Activity Projection Difference";"BMI-Ball"])
xticks([1,2])
xticklabels(["Raw W";"Absolute W"])

subplot(1,2,2)
yline(0,'--','LineWidth',2);
title("View Angle Velocity")
ylabel(["Mean Activity Projection Difference";"BMI-Ball"])
xticks([1,2])
xticklabels(["Raw W";"Absolute W"])

figure
i = 1;
scatter(plot_off+(i)*ones(1,num_sess),(res_mat_on(i,2,:)-res_mat_on(i,1,:)),'filled','k')
hold on
plot([plot_off(1)+(i),plot_off(end)+(2)],[mean((res_mat_on(i,2,:)-res_mat_on(i,1,:)),'omitnan'),mean((res_mat_on(i,2,:)-res_mat_on(i,1,:)),'omitnan')],'k','LineWidth',2)  

scatter(plot_off+(2)*ones(1,num_sess),(res_mat_off(i,2,:)-res_mat_off(i,1,:)),'filled','k')
hold on
plot([plot_off(1)+(2),plot_off(end)+(2)],[mean((res_mat_off(i,2,:)-res_mat_off(i,1,:)),'omitnan'),mean((res_mat_off(i,2,:)-res_mat_off(i,1,:)),'omitnan')],'k','LineWidth',2) 

yline(0,'--','LineWidth',2);
% title("View Angle Velocity")
ylabel(["Mean Activity Projection Difference";"BMI-Ball"])
xticks([1,2])
xticklabels(["View Angle";"View Angle Velocity"])
xtickangle(30)

%% Plot scatter va va vav 

figure

scatter((res_mat_on(1,2,:)-res_mat_on(1,1,:))./res_mat_on(1,1,:),(res_mat_off(1,2,:)-res_mat_off(1,1,:))./res_mat_off(1,1,:),'filled','k')
hold on

min_val_on = min((res_mat_on(1,2,:)-res_mat_on(1,1,:))./res_mat_on(1,1,:),[],'all');
min_val_off = min((res_mat_off(1,2,:)-res_mat_off(1,1,:))./res_mat_off(1,1,:),[],'all');
max_val_on = max((res_mat_on(1,2,:)-res_mat_on(1,1,:))./res_mat_on(1,1,:),[],'all');
max_val_off = max((res_mat_off(1,2,:)-res_mat_off(1,1,:))./res_mat_off(1,1,:),[],'all');

min_val = min([min_val_on,min_val_off]);
max_val = min([max_val_on,max_val_off]);

% plot([min_val,max_val],[min_val,max_val],'--','LineWidth',2,'Color',[0.5,0.5,0.5])
plot([-0.4,0.4],[-0.4,0.4],'--','LineWidth',2,'Color',[0.5,0.5,0.5])

xlabel("View Angle")
ylabel("View Angle Velocity")
title(["Mean Activity Projection Difference";"BMI-Ball"])

xline(0,'--','LineWidth',2);
yline(0,'--','LineWidth',2);
axis 'square'
xlim([-0.4,0.4])
ylim([-0.4,0.4])

%% alternative z-scoring

figure

scatter((res_mat_on(1,2,:)-res_mat_on(1,1,:))./std_mat_on(1,1,:),(res_mat_off(1,2,:)-res_mat_off(1,1,:))./std_mat_off(1,1,:),'filled','k')
hold on

min_val_on = min((res_mat_on(1,2,:)-res_mat_on(1,1,:))./std_mat_on(1,1,:),[],'all');
min_val_off = min((res_mat_off(1,2,:)-res_mat_off(1,1,:))./std_mat_off(1,1,:),[],'all');
max_val_on = max((res_mat_on(1,2,:)-res_mat_on(1,1,:))./std_mat_on(1,1,:),[],'all');
max_val_off = max((res_mat_off(1,2,:)-res_mat_off(1,1,:))./std_mat_off(1,1,:),[],'all');

min_val = min([min_val_on,min_val_off]);
max_val = min([max_val_on,max_val_off]);

% plot([min_val,max_val],[min_val,max_val],'--','LineWidth',2,'Color',[0.5,0.5,0.5])
plot([-0.5,0.5],[-0.5,0.5],'--','LineWidth',2,'Color',[0.5,0.5,0.5])

xlabel("View Angle")
ylabel("View Angle Velocity")
title(["Mean Activity Projection Difference";"BMI-Ball"])

xline(0,'--','LineWidth',2);
yline(0,'--','LineWidth',2);
axis 'square'
xlim([-0.5,0.5])
ylim([-0.5,0.5])

%% Heirarchical bootstrap
boot_samps = 1000;
num_trials = 4;
% Calculates probability that data2 >= data1. 

%Make nan final value in row
res_mat_on(:,:,4,[2,3]) = res_mat_on(:,:,4,[3,4]);
res_mat_on(:,:,4,4) = nan;

res_mat_off(:,:,4,[2,3]) = res_mat_off(:,:,4,[3,4]);
res_mat_off(:,:,4,4) = nan;
% 
% [p_boot_1, bootstats, bootstats_center, bootstats_sem] = get_bootstrap_results_equalsamples(squeeze((res_mat_on(1,2,:,:)-res_mat_on(1,1,:,:))./res_mat_on(1,1,:,:)),squeeze((res_mat_off(1,2,:,:)-res_mat_off(1,1,:,:))./res_mat_off(1,1,:,:)),boot_samps,num_trials,'mean');
% [p_boot_2, bootstats, ~, ~] = get_bootstrap_results_equalsamples(squeeze((res_mat_off(1,2,:,:)-res_mat_off(1,1,:,:))./res_mat_off(1,1,:,:)),squeeze((res_mat_on(1,2,:,:)-res_mat_on(1,1,:,:))./res_mat_on(1,1,:,:)),boot_samps,num_trials,'mean');

% z-scored
boot_samps = 1000;
num_trials = 4;
% Calculates probability that data2 >= data1. 

%Make nan final value in row - not best way
std_mat_on(:,:,4,[2,3]) = std_mat_on(:,:,4,[3,4]);
std_mat_on(:,:,4,4) = nan;

std_mat_off(:,:,4,[2,3]) = std_mat_off(:,:,4,[3,4]);
std_mat_off(:,:,4,4) = nan;

% [p_boot_1, bootstats, bootstats_center, bootstats_sem] = get_bootstrap_results_equalsamples(squeeze((res_mat_on(1,2,:,:)-res_mat_on(1,1,:,:))./std_mat_on(1,1,:,:)),squeeze((res_mat_off(1,2,:,:)-res_mat_off(1,1,:,:))./std_mat_off(1,1,:,:)),boot_samps,num_trials,'mean');
% [p_boot_2, bootstats, ~, ~] = get_bootstrap_results_equalsamples(squeeze((res_mat_off(1,2,:,:)-res_mat_off(1,1,:,:))./std_mat_off(1,1,:,:)),squeeze((res_mat_on(1,2,:,:)-res_mat_on(1,1,:,:))./std_mat_on(1,1,:,:)),boot_samps,num_trials,'mean');

[all_p_boot,all_centres,all_sems] = run_H_boot_ets(squeeze((res_mat_on(1,2,:,:)-res_mat_on(1,1,:,:))./std_mat_on(1,1,:,:)), squeeze((res_mat_off(1,2,:,:)-res_mat_off(1,1,:,:))./std_mat_off(1,1,:,:)),false);

h_boots.all_p_boot = all_p_boot;
h_boots.all_centres = all_centres;
h_boots.all_sems = all_sems;