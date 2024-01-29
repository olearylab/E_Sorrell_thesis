function [wilc_res,h_boots] = plot_pop_mean_results_06122022(all_popboot_mean,all_overall_n_means,centres)
% 06/12/2022

% plot results for mean value comparisons across trial types.

% edits for paper ready figures

num_mice = size(all_popboot_mean,1);
num_days = size(all_popboot_mean,2);
nbins = length(centres);
% colour_vec = [[0 0.4470 0.7410];[0.8500 0.3250 0.0980];[0.9290 0.6940 0.1250];[0.4940 0.1840 0.5560];[0.4660 0.6740 0.1880];[0.3010 0.7450 0.9330];[0.6350 0.0780 0.1840]];
colour_vec = [[0 0 0];[0 0 0];[0 0 0];[0 0 0];[0 0 0];[0 0 0];[0 0 0]];

% centres = centres*0.74;
%% Plot mean changes by bin

% Plot average across mice and sessions, with faint lines for average for
% individual mice (could change to faint line for each session).

% Left and right separately
% figure
% mean of differences for each mouse
% m_means = nan.*ones(num_mice,2,length(centres));
% % ready for H bootstrapping
% md_means = nan.*ones(num_mice,num_days,2,length(centres));
% for m = 1:num_mice
%     cur_m = nan.*ones(num_days,2,length(centres));
%     for d = 1:num_days
%         if ~isempty(all_popboot_mean{m,d})
%             cur_means = all_popboot_mean{m,d};
%             for i = 1:2
%                 cur_m(d,i,:) = cur_means(i+2,:)-cur_means(i,:);
%             end           
%             
%         end
%         
%     end
%     m_means(m,:,:) = squeeze(mean(cur_m,'omitnan'));
%     md_means(m,1:length(~isnan(cur_m(:,1,1))),:,:) = cur_m(~isnan(cur_m(:,1,1)),:,:);
% end

% for i = 1:2
%     subplot(1,2,i)
%     yline(0,'LineWidth',2);
%     hold on
% end

% for m = 1:num_mice
%     for i = 1:2
%         subplot(1,2,i)
%         plot(centres,squeeze(m_means(m,i,:)),'-','LineWidth',1.5,'Color',colour_vec(m,:))
%         hold on
%     end
% end
% titles = ["Left";"Right"];
% for i = 1:2
%     subplot(1,2,i)
%     plot(centres,squeeze(mean(m_means(:,i,:),'omitnan')),'LineWidth',3,'Color','k')
%     xlabel("Linearised Position (cm)")
%     xticks([0,200])
%     ylim([-0.4,0.6])
%     title(titles(i))
%     box off
% end

% subplot(1,2,1)
% ylabel(["Mean Activity Difference (a.u.)";"BMI - Ball"])
% Average for both directions:
figure
% subplot(1,3,3)
% mean of differences for each mouse
m_means = nan.*ones(num_mice,length(centres));

md_means = nan.*ones(num_mice,num_days,length(centres));
for m = 1:num_mice
    cur_m = nan.*ones(num_days,length(centres));
    for d = 1:num_days
        if ~isempty(all_popboot_mean{m,d})
            cur_means = all_popboot_mean{m,d};
            av_lr = zeros(2,length(centres));
            for i = 1:2
                av_lr(i,:) = mean(cur_means([(i-1)*2+1,i*2],:));
            end  
            cur_m(d,:) = av_lr(2,:)- av_lr(1,:);
          
            
        end
        
    end
    m_means(m,:) = squeeze(mean(cur_m,'omitnan'));
    md_means(m,1:sum(~isnan(cur_m(:,1))),:) = cur_m(~isnan(cur_m(:,1)),:);
end

yline(0,'--','LineWidth',2);
hold on

for m = 1:num_mice

    plot(centres,squeeze(m_means(m,:)),'-','LineWidth',1.5,'Color',[0.5,0.5,0.5])

end


plot(centres,squeeze(mean(m_means(:,:),'omitnan')),'LineWidth',3,'Color','k')
xlabel("Linearised Position (cm)")
xticks([0,200])
ylabel(["Mean Normalised Activity";"Differece, BMI - Ball (a.u.)"])
box off
axis('square')

figure
% bootstrapped mean and std
boot_samps = 1000;
num_trials = 4;
% Calculates probability that data2 > data1.
bootstats = NaN(nbins,boot_samps);
for b = 1:nbins
    bootstats(b,:) = get_bootstrapped_equalsamples(squeeze(md_means(:,:,b)),boot_samps,num_trials,'mean');
end

% CI_vals = [0.5,99.5]; % [2.5,97.5];
CI_vals = [2.5,97.5];

% CIs, upper/lower x ntypes x nbins x zdim
CIs_all = zeros(2,nbins);
for i = 1:nbins
    CIs_all(:,i) = prctile(bootstats(i,:),CI_vals);
end

% boostrap stds instead
boot_stds = std(bootstats,0,2);
lims_all = zeros(2,nbins);
lims_all(1,:) = mean(bootstats,2,'omitnan') - boot_stds;
lims_all(2,:) = mean(bootstats,2,'omitnan') + boot_stds;

% h = fill([centres,fliplr(centres)],[CIs_all(1,:),fliplr(CIs_all(2,:))],'k','EdgeColor','none');
h = fill([centres,fliplr(centres)],[lims_all(1,:),fliplr(lims_all(2,:))],'k','EdgeColor','none');
set(h,'facealpha',.3)
hold on

plot(centres,squeeze(mean(bootstats,2,'omitnan')),'LineWidth',2,'Color','k')

xlabel("Linearised Position (cm)")
xticks([0,200])
ylabel(["Mean Normalised Activity";"Differece, BMI - Ball (a.u.)"])
box off
axis('square')
yline(0,'--','LineWidth',2);
xline(200*0.74,'--','Linewidth',2); % Cue end
xline(300*0.74,'--','Linewidth',2); % Turn start

%% Plot overall mean changes by neuron

n_diff_cell = cell(num_mice,1);
% calculate cumulative distribution for each mouse (could do each session
% separately, then average acros sessions? Average proportions...
num_points = 1000;
cum_line = zeros(num_mice,num_points);
m_vals = zeros(num_mice,num_points);

% prepare for h bootstrapping. proportion increasing
prop_ready = nan.*ones(num_mice,num_days);
for m = 1:num_mice
    cur_m_diffs = [];
    d_ind = 0;
    for d = 1:num_days
        if ~isempty(all_overall_n_means{m,d})
            d_ind = d_ind+1;
            cur_means = all_overall_n_means{m,d};

            cur_diffs = cur_means(2,:)-cur_means(1,:);
            
            cur_m_diffs = [cur_m_diffs,cur_diffs];
            
            prop_ready(m,d_ind) = sum(cur_diffs>0)/length(cur_diffs);
            
        end
    end
    n_diff_cell{m} = cur_m_diffs;
    cur_vals = linspace(min(cur_m_diffs),max(cur_m_diffs),num_points);
    m_vals(m,:) = cur_vals;
    for i = 1:num_points
        cum_line(m,i) = sum(cur_m_diffs<cur_vals(i))/length(cur_m_diffs);
    end
end
% figure
% violin(n_diff_cell);
% yline(0,'LineWidth',2);

% bars showing change

all_diffs = [];
for m = 1:num_mice
    for d = 1:num_days
        if ~isempty(all_overall_n_means{m,d})
            cur_means = all_overall_n_means{m,d};

            cur_diffs = cur_means(2,:)-cur_means(1,:);
            
            all_diffs = [cur_m_diffs,cur_diffs];
            
        end
    end
end
all_vals = linspace(min(all_diffs),max(all_diffs),num_points);
all_cum_line = zeros(1,num_points);
for i = 1:num_points
    all_cum_line(i) = sum(all_diffs<all_vals(i))/length(all_diffs);
end

rng(1);
boot_samps = 1000;
num_trials = 4;
bootstats = get_bootstrapped_equalsamples(prop_ready,boot_samps,num_trials,'mean');
%Get mean and SEM of bootstrapped samples:
h_boots.all_sems = std(bootstats);
h_boots.all_centres = mean(bootstats);
% figure
% barh(sort(all_diffs))

% cumulative line
% figure
% 
% cdfplot(all_diffs)
% xline(0,'LineWidth',2);

% cumulative line - self constructed
figure
% subplot(1,3,2)

for m = 1:num_mice
    plot(m_vals(m,:),cum_line(m,:),'-','LineWidth',1.5,'Color',[0.5,0.5,0.5])
    hold on
end

plot(all_vals,all_cum_line,'-','LineWidth',3,'Color','k')
xline(0,'--','LineWidth',2);
yline(0.5,'--','LineWidth',2);
box off
axis('square')
xlabel(["Mean Normalised Activity";"Difference, BMI - Ball (a.u.)"])
ylabel(["Cumulative Fraction"; "of Neurons"])

% figure
% violin(all_diffs');

% figure
% histogram(all_diffs);
% xline(0,'LineWidth',2);

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

%% Calculate significance of difference with two-sided Wilcoxon signed rank
wilc_res = signrank(full_means(1,:),full_means(2,:));

%% Plot overall means - alternative figure ready version averaged

% Edited for single bar

figure
% subplot(1,3,1)
bar(overall_means(2)-overall_means(1),'w','LineWidth',2)
hold on

for m = 1:num_mice
    plot((m_means(2,m)-m_means(1,m)),'o','LineWidth',1.5,'Color',colour_vec(m,:),'MarkerFaceColor',colour_vec(m,:))
end
xticks(1)
xticklabels(["BMI - Ball"])
% xlabel("Trial Type")
ylim([0,0.4])
yticks([0,0.2,0.4])
ylabel(["Mean Normalised Activity";"Difference (a.u.)"])
box off
axis('square')

%% Plot overall means over time
% figure
% subplot(1,2,1)
% for m = 1:num_mice
%     plot(squeeze(full_means(1,m,:)),'-o','LineWidth',1.5,'Color',colour_vec(m,:),'MarkerFaceColor',colour_vec(m,:))
%     hold on
%     plot(squeeze(full_means(2,m,:)),'--o','LineWidth',1.5,'Color',colour_vec(m,:),'MarkerFaceColor',colour_vec(m,:))
% end
% % xticks([1,2])
% % xticklabels(["Ball";"BMI"])
% xlabel("Day")
% % ylim([-0.25,0.25])
% % yticks([-0.25,0,0.25])
% ylabel("Mean Normalised Activity (a.u.)")
% box off
% 
% subplot(1,2,2)
% for m = 1:num_mice
%     plot(squeeze(full_means(2,m,:))-squeeze(full_means(1,m,:)),'-o','LineWidth',1.5,'Color',colour_vec(m,:),'MarkerFaceColor',colour_vec(m,:))
%     hold on
% end
% % xticks([1,2])
% % xticklabels(["Ball";"BMI"])
% xlabel("Day")
% % ylim([-0.25,0.25])
% % yticks([-0.25,0,0.25])
% ylabel(["Mean Activity Difference"; "BMI - Ball (a.u.)"])
% box off