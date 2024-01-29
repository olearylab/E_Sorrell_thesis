function [] = plot_across_sessions_mean_change_18042023(all_overall_n_means)
% 18/04/2023

% Plot mean changes by day to see if there is a trend over time. Does it
% matter than each day is normalised separately? Possibly. Could re-run
% other code to normalise all together if desired.

num_mice = size(all_overall_n_means,1);
num_days = size(all_overall_n_means,2);
% colour_vec = [[0 0.4470 0.7410];[0.8500 0.3250 0.0980];[0.9290 0.6940 0.1250];[0.4940 0.1840 0.5560];[0.4660 0.6740 0.1880];[0.3010 0.7450 0.9330];[0.6350 0.0780 0.1840]];
colour_vec = [[0 0 0];[0 0 0];[0 0 0];[0 0 0];[0 0 0];[0 0 0];[0 0 0]];

%% Plot overall means
full_means = nan.*ones(2,num_mice,num_days);
for m = 1:num_mice
    for d = 1:num_days
        if ~isempty(all_overall_n_means{m,d})
            cur_means = all_overall_n_means{m,d};
            cur_means = mean(cur_means,2,'omitnan');
            full_means(:,m,d) = cur_means;         
            
        end
    end
end
overall_means = squeeze(mean(full_means,[2,3],'omitnan'));



m_means = squeeze(mean(full_means,3,'omitnan'));

%% Plot overall means over time
figure

for m = 1:num_mice
    cur_res = squeeze(full_means(2,m,:))-squeeze(full_means(1,m,:));
    plot(find(~isnan(cur_res)),cur_res(~isnan(cur_res)),'-o','LineWidth',1.5,'Color',[0.5,0.5,0.5],'MarkerFaceColor',[0.5,0.5,0.5])
    hold on
end
plot(mean(squeeze(full_means(2,:,:))-squeeze(full_means(1,:,:)),1,'omitnan'),'-o','LineWidth',3,'Color','k','MarkerFaceColor','k')
% xticks([1,2])
% xticklabels(["Ball";"BMI"])
xlabel("Day")
% ylim([-0.25,0.25])
% yticks([-0.25,0,0.25])
ylabel(["Mean Activity Difference"; "BMI - Ball (a.u.)"])
title(["Mean Activity Difference";"Across Sessions"])
box off
axis('square')
xticks([1,2,3,4])
xlim([0.5,4.5])