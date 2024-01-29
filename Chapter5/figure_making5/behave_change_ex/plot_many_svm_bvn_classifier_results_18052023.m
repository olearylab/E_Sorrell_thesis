function [] = plot_many_svm_bvn_classifier_results_18052023(accuracies_cell_bvn,centres)
% 17/02/2023

% Function for plotting results of svm classification of trials as left or
% right using behaviour (pitch and yaw)

num_mice = size(accuracies_cell_bvn,1);
num_days = size(accuracies_cell_bvn,2);

means_all = nan.*ones(num_mice,num_days,2,4,length(centres));

centres = centres*0.74;

% prepare for summary figure
for m = 1:num_mice
    for d = 1:num_days
        if ~isempty(accuracies_cell_bvn{m,d})
            cur_acc = accuracies_cell_bvn{m,d};
            
            for i = 1:2
                means_all(m,d,i,1,:) = cur_acc(i,:);
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

%% Plot combined
colour_vec = [[0.4940 0.1840 0.5560];[0.4660 0.6740 0.1880]];
titles = ["Ball vs BMI Classifier"];
figure
cur_mean = permute(squeeze(mean(squeeze(means_all(:,:,:,1,:)),3,'omitnan')),[3,1,2]);
plot(centres,cur_mean(:,:),'Color',[0.5,0.5,0.5])
hold on
% errorbar(centres,mean(cur_mean,[2,3],'omitnan'),std(cur_mean,0,[2,3],'omitnan'),'LineWidth',2,'Color',colour_vec(j))
plot(centres,mean(cur_mean,[2,3],'omitnan'),'LineWidth',3,'Color','k')
yline(0.5,'--','Linewidth',2);
xline(200*0.74,'--','Linewidth',2); % Cue end
xline(300*0.74,'--','Linewidth',2); % Turn start
title(titles)
ylim([0,1])
xlabel("Linearised Position (cm)")

ylabel(["Classification Accuracy"])
axis('square')