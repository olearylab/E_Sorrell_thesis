function [h_boots] = plot_many_svm_bvn_neural_classifier_results(accuracies_cell_bvn,accuracies_cell_shuff,centres)
% 03/08/2023

% Function for plotting results of svm classification of trials as ball vs 
% bmi using neural data

% Will need to reassess what error bars to use etc.

num_mice = size(accuracies_cell_bvn,1);
num_days = size(accuracies_cell_bvn,2);

means_all = nan.*ones(num_mice,num_days,2,4,length(centres));
means_all_shuff = nan.*ones(num_mice,num_days,2,4,length(centres));

centres = centres*0.74;

nbins = length(centres);

% prepare for summary figure
for m = 1:num_mice
    for d = 1:num_days
        if ~isempty(accuracies_cell_bvn{m,d})
            cur_acc = accuracies_cell_bvn{m,d};
            cur_acc_shuff = accuracies_cell_shuff{m,d};
            
            for i = 1:2
                means_all(m,d,i,1,:) = cur_acc(i,:);
                means_all_shuff(m,d,i,1,:) = mean(cur_acc_shuff(:,i,:));
            end

        end
    end
end

% Plot
% % colour_vec = [[0.4940 0.1840 0.5560];[0.4660 0.6740 0.1880]];
% % titles = ["ball vs bmi classifier"];
% % figure
% % for j = 1:2
% %     cur_mean = permute(squeeze(means_all(:,:,j,1,:)),[3,1,2]);
% %     plot(centres,cur_mean(:,:),'Color',colour_vec(j,:))
% %     hold on
% %     % errorbar(centres,mean(cur_mean,[2,3],'omitnan'),std(cur_mean,0,[2,3],'omitnan'),'LineWidth',2,'Color',colour_vec(j))
% %     plot(centres,mean(cur_mean,[2,3],'omitnan'),'LineWidth',3,'Color',colour_vec(j,:))
% % 
% % end
% % yline(0.5,'--','Linewidth',2);
% % xline(200*0.74,'--','Linewidth',2); % Cue end
% % xline(300*0.74,'--','Linewidth',2); % Turn start
% % title(titles)
% % ylim([0,1])
% % xlabel("Linearised Position (cm)")
% % 
% % ylabel(["Classification Accuracy"])
% % axis('square')
% % 
% % % Plot - with std across sessions rather than all session lines.
% % % Actually use 95% confidence interval rather than std
% % colour_vec = [[0.4940 0.1840 0.5560];[0.4660 0.6740 0.1880]];
% % titles = ["ball vs bmi classifier"];
% % figure
% % for j = 1:2
% %     cur_mean = permute(squeeze(means_all(:,:,j,1,:)),[3,1,2]);
% % 
% %     CIs = prctile(cur_mean(:,:),CI_vals,2);
% % 
% %     plot(centres,CIs(:,1),'Linewidth',1,'Color',colour_vec(j,:))
% %     hold on
% %     plot(centres,CIs(:,2),'Linewidth',1,'Color',colour_vec(j,:))
% % 
% %     h = fill([centres,fliplr(centres)],[CIs(:,1)',fliplr(CIs(:,2)')],colour_vec(j,:));
% %     set(h,'facealpha',.1)
% % 
% %     % errorbar(centres,mean(cur_mean,[2,3],'omitnan'),std(cur_mean,0,[2,3],'omitnan'),'LineWidth',2,'Color',colour_vec(j))
% %     plot(centres,mean(cur_mean,[2,3],'omitnan'),'LineWidth',3,'Color',colour_vec(j,:))
% %     
% % 
% % end
% % 
% % yline(0.5,'--','Linewidth',2);
% % xline(200*0.74,'--','Linewidth',2); % Cue end
% % xline(300*0.74,'--','Linewidth',2); % Turn start
% % title(titles)
% % ylim([0,1])
% % xlabel("Linearised Position (cm)")
% % 
% % ylabel(["Classification Accuracy"])
% % axis('square')

num_sess = 19;
%% Plot combined
colour_vec = [[0.4940 0.1840 0.5560];[0.4660 0.6740 0.1880]];
% titles = ["Ball vs BMI Classifier"];
figure
cur_mean = permute(squeeze(mean(squeeze(means_all(:,:,:,1,:)),3,'omitnan')),[3,1,2]);
cur_mean_shuff = permute(squeeze(mean(squeeze(means_all_shuff(:,:,:,1,:)),3,'omitnan')),[3,1,2]);
% 
% cur_sem = std(cur_mean,0,[2,3],'omitnan');% ./sqrt(num_sess);
% cur_sem_shuff = std(cur_mean_shuff,0,[2,3],'omitnan');% ./sqrt(num_sess);

% H bootstrap means and stds
boot_samps = 1000;
num_trials = 4;
% Calculates probability that data2 > data1. 
% [p_boot, bootstats, bootstats_center, bootstats_sem] = get_bootstrap_results_equalsamples(squeeze(p_R2_all_mat(:,:,2)),squeeze(p_R2_all_mat(:,:,4)),boot_samps,num_trials,'mean');
orig_cur_means = cur_mean;
for m = 1:num_mice
    d_ind = 0;
    for d = 1:num_days
        if ~isnan(orig_cur_means(1,m,d)) 
            d_ind = d_ind+1;
            cur_mean(:,m,d_ind) = orig_cur_means(:,m,d);

        end
    end
    if d_ind<num_days
        cur_mean(:,m,d_ind+1:num_days) = nan;
    end
end

orig_cur_mean_shuff = cur_mean_shuff;
for m = 1:num_mice
    d_ind = 0;
    for d = 1:num_days
        if ~isnan(orig_cur_mean_shuff(1,m,d)) 
            d_ind = d_ind+1;
            cur_mean_shuff(:,m,d_ind) = orig_cur_mean_shuff(:,m,d);

        end
    end
    if d_ind<num_days
        cur_mean_shuff(:,m,d_ind+1:num_days) = nan;
    end
end

all_centres = NaN(nbins,2);
all_sems = NaN(nbins,2);
all_p_boot = NaN(nbins,1);
for b = 1:nbins
    [all_p_boot(b),all_centres(b,:),all_sems(b,:)] = run_H_boot_ets(squeeze(cur_mean(b,:,:)), squeeze(cur_mean_shuff(b,:,:)),false);
end

h_boots.all_p_boot = all_p_boot;
h_boots.all_centres = all_centres;
h_boots.all_sems = all_sems;

% bootstats_s = NaN(nbins,boot_samps);
% for b = 1:nbins
%     bootstats_s(b,:) = get_bootstrapped_equalsamples(squeeze(cur_mean_shuff(b,:,:)),boot_samps,num_trials,'mean');
% end
% 
% bootstats = NaN(nbins,boot_samps);
% for b = 1:nbins
%     bootstats(b,:) = get_bootstrapped_equalsamples(squeeze(cur_mean(b,:,:)),boot_samps,num_trials,'mean');
% end

% CI_vals = [0.5,99.5]; % [2.5,97.5];
% CI_vals = [2.5,97.5];
% 
% % CIs, upper/lower x ntypes x nbins x zdim
% CIs_all = zeros(2,nbins);
% for i = 1:nbins
%     CIs_all(:,i) = prctile(bootstats(i,:),CI_vals);
% end

% boostrap stds instead
boot_stds = all_sems;
lims_all = zeros(2,nbins,2);
lims_all(1,:,:) = all_centres - boot_stds;
lims_all(2,:,:) = all_centres + boot_stds;


% plot(centres,cur_mean(:,:),'Color',[0.5,0.5,0.5])

% h = fill([centres,fliplr(centres)],[CIs_all(1,:),fliplr(CIs_all(2,:))],'k','EdgeColor','none');
h = fill([centres,fliplr(centres)],[squeeze(lims_all(1,:,1)),fliplr(squeeze(lims_all(2,:,1)))],'k','EdgeColor','none');
set(h,'facealpha',.3)
hold on

plot(centres,all_centres(:,1),'LineWidth',2,'Color','k')

% CI_vals = [0.5,99.5]; % [2.5,97.5];
% CI_vals = [2.5,97.5];
% 
% % CIs, upper/lower x ntypes x nbins x zdim
% CIs_all_s = zeros(2,nbins);
% for i = 1:nbins
%     CIs_all_s(:,i) = prctile(bootstats_s(i,:),CI_vals);
% end

% boostrap stds instead
% boot_stds_s = std(bootstats_s,0,2);
% lims_all_s = zeros(2,nbins);
% lims_all_s(1,:) = mean(bootstats_s,2,'omitnan') - boot_stds_s;
% lims_all_s(2,:) = mean(bootstats_s,2,'omitnan') + boot_stds_s;

% h = fill([centres,fliplr(centres)],[CIs_all_s(1,:),fliplr(CIs_all_s(2,:))],'k','EdgeColor','none');
h = fill([centres,fliplr(centres)],[squeeze(lims_all(1,:,2)),fliplr(squeeze(lims_all(2,:,2)))],'k','EdgeColor','none');
set(h,'facealpha',.3)
hold on

plot(centres,all_centres(:,2),'--','LineWidth',2,'Color','k')


% yline(0.5,'--','Linewidth',2);
xline(200*0.74,'--','Linewidth',2); % Cue end
xline(300*0.74,'--','Linewidth',2); % Turn start
title(["Ball vs BMI Classifier";"Using Neural Data"])
ylim([0,1])
xlabel("Linearised Position (cm)")

ylabel(["Classification Accuracy"])
axis('square')

%% get p_boots
% all_p_boot = nan.*ones(nbins,1);
% for b = 1:nbins
%     all_p_boot(b) = get_direct_prob(bootstats_s(b,:),bootstats(b,:));
% end
    