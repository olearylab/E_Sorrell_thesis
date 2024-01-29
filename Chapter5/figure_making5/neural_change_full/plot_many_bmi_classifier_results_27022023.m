function [] = plot_many_bmi_classifier_results_27022023(x_binned_all,x_binned_all_r,virmen_cell)
% 27/02/2023

% Function for plotting results of bmi vs normal classifier for many mice
% and days.

% Take results for left trials and for right trials. Create separate left
% and right figure, and create overall combined figure.

% Need to rethink calculating the average across left/right (once these are
% properly balanced it should be fine)

num_mice = size(x_binned_all,1);
num_days = size(x_binned_all,2);
nbins = size(x_binned_all{1,1},2);

virmen_ex = virmen_cell{1,1};
virmen_ex(6,:) = virmen_ex(6,:) + abs(virmen_ex(5,:));
[~,edges] = discretize(virmen_ex(6,virmen_ex(8,:)==0),nbins);

%% Plot results averaged across days for each trial type 
% averaging across days is false average as not scaled by trial type. Need
% to be more careful about all this. Maybe save number of trials, or work
% with all results. Maybe go back to original method.
% t_types = [1,2,3,4];
% types_vec = [1,4,7,10];
% types_vec = types_vec(t_types);
final_results = zeros(num_mice,2,nbins);
final_results_r = zeros(num_mice,2,nbins);
for m = 1:num_mice
    m_binned = nan.*ones(num_days,2,nbins);
    m_binned_r = nan.*ones(num_days,2,nbins);
    for d = 1:num_days
        if ~isempty(x_binned_all{m,d})
            m_binned(d,:,:) = x_binned_all{m,d};
            m_binned_r(d,:,:)= x_binned_all_r{m,d};
        end
    end
    final_results(m,:,:) = mean(m_binned,'omitnan');
    final_results_r(m,:,:) = mean(m_binned_r,'omitnan');
end

lr_results = zeros(size(final_results));

lr_results(:,1,:) = squeeze(mean(final_results,2));
lr_results(:,2,:) = squeeze(mean(final_results_r,2));


%% Plotting
% Left and right separate
pos_scale = 0.74;
centres = pos_scale*(edges(2:end)+edges(1:end-1))/2;
edges = edges*pos_scale;
cue_end = 200*pos_scale;
turn_point = 300*pos_scale; % Could pick 305   

% plot_types_c = [1,4,7,10];
% % titles_c = ["Left";"Right"];
% % 
% % figure
% % 
% % for j = 1:2
% %     subplot(2,1,j)
% %     hold on
% %     xline(turn_point,'LineWidth',2,'Color','k');
% %     xline(cue_end,'LineWidth',2,'Color','k');
% %     yline(0.5,'LineWidth',2);
% %     
% %     for m = 1:num_mice
% %         plot(centres,squeeze(lr_results(m,j,:)),'LineWidth',2)
% %     end
% %     plot(centres,squeeze(mean(lr_results(:,j,:))),'LineWidth',2,'Color','k')
% % 
% %     title(titles_c(j))
% %     ylim([0,1])
% %     if (j == 1) || (j == 3)
% %         ylabel(["Classification"; "Accuracy"])
% %         yticks([0,0.5,1])
% %     else
% %         yticks([])
% %     end
% %     if j > 2
% %         xlabel("Position (cm)")
% %         xticks([0,200])
% %     else
% %         xticks([])
% %     end
% % end

% Plot average over left and right

av_results = squeeze(mean(lr_results,2));

figure

hold on
xline(turn_point,'--','LineWidth',2,'Color','k');
xline(cue_end,'--','LineWidth',2,'Color','k');
yline(0.5,'--','LineWidth',2);
    
for m = 1:num_mice
    plot(centres,squeeze(av_results(m,:)),'LineWidth',2,'Color',[0.5,0.5,0.5])
end
plot(centres,squeeze(mean(av_results)),'LineWidth',3,'Color','k')

title(["BMI vs Ball"; "Neural Classification"])
ylim([0,1])

ylabel(["Classification"; "Accuracy"])
yticks([0,0.5,1])

xlabel("Linearised Position (cm)")
axis('square')
