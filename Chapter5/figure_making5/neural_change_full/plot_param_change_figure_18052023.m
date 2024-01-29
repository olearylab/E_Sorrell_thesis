function [h_boots] = plot_param_change_figure_18052023(sig_results_all,sig_summary_cell,all_full_results_both,full_all_param_means,sig_check_all,bootstrap_means,centres,ex_neu,ex_dir,ex_md)
% 18/05/2023

% Need to pick suitable colours
% colours_vec = [[0 0.4470 0.7410];[0.6350 0.0780 0.1840];[0.9290 0.6940 0.1250];[0.4660 0.6740 0.1880]];
% colours_vec = [[0 0.4470 0.7410];[0.8500 0.3250 0.0980]];
colours_vec = [[0.4940 0.1840 0.5560];[0.4660 0.6740 0.1880]];
centres = centres*0.74;
CI_vals = [2.5,97.5];
t_types = [1,4,7,10];
%% Now summary plots
figure

num_mice = size(all_full_results_both,1);
num_days = size(all_full_results_both,2);

% Want proportions of: gain, loss, both no change, pos, amp, width
all_proportions_both = nan.*ones(num_mice,num_days,2,6);
% all_proportions_both_dir = nan.*ones(num_mice,num_days,6);

num_tuned = zeros(num_mice,num_days,4);

for m = 1:num_mice
    for d = 1:num_days
        if ~isempty(all_full_results_both{m,d})
            full_results_both = all_full_results_both{m,d};
            cur_sig = sig_results_all{m,d};
            cur_sig_summary = sig_summary_cell{m,d};
            total_tuned = sum(cur_sig==1)>0;
            for i = 1:2
                cur_neurons = sum(cur_sig([i,i+2],:)==1)>0;
                cur_both = squeeze(cur_sig_summary(i,1,:)==1);
                all_proportions_both(m,d,i,1) = sum(squeeze(cur_sig_summary(i,2,:))==1)./sum(cur_neurons);
                all_proportions_both(m,d,i,2) = sum(squeeze(cur_sig_summary(i,3,:))==1)./sum(cur_neurons);
                % Change to sum(cur_neurons) instead of sum(cur_both) to get
                % proportion of tuned to left/right matching above.
                all_proportions_both(m,d,i,3) = squeeze(sum(sum(squeeze(full_results_both(cur_both,i,:))==0,2)==3))./sum(cur_neurons);
                all_proportions_both(m,d,i,4:6) = squeeze(sum(squeeze(full_results_both(cur_both,i,:))))./sum(cur_neurons);
            end
            % average over left and right, rather than total with at least
            % one and dividing by tuned to at least one
            % cur_neurons = sum(cur_sig==1)>0;

            % Alternatively, get proportion of total neurons tuned to
            % something that fit each description
            % all_proportions_both_dir(m,d,1) = sum(squeeze(sum(cur_sig_summary(:,2,total_tuned))>0))./sum(total_tuned);
            % all_proportions_both_dir(m,d,2) = sum(squeeze(sum(cur_sig_summary(:,3,total_tuned))>0))./sum(total_tuned);
%             
            % cur_both_l = squeeze(cur_sig_summary(1,1,:)==1);
            % cur_both_r = squeeze(cur_sig_summary(2,1,:)==1);
%             % No change in either direction
            % all_proportions_both_dir(m,d,3) = squeeze(sum(sum(sum(full_results_both(cur_both_l & cur_both_r,:,:)==0,3)==3,2)==2))./sum(total_tuned);
            % No change in at least on direction
            % all_proportions_both_dir(m,d,3) = (squeeze(sum(sum(full_results_both(cur_both_l,1,:)==0,3)==3,2)) + squeeze(sum(sum(full_results_both(cur_both_r,2,:)==0,3)==3,2)) - squeeze(sum(sum(sum(full_results_both(cur_both_l & cur_both_r,:,:)==0,3)==3,2)==2)))./sum(total_tuned);
            
            % all_proportions_both_dir(m,d,4:6) = (squeeze(sum(full_results_both(cur_both_l,1,:))) + squeeze(sum(full_results_both(cur_both_r,2,:))) - squeeze(sum(sum(full_results_both(cur_both_l & cur_both_r,:,:),2)==2)))./sum(total_tuned);
            
            for i = 1:4
                num_tuned(m,d,i) = sum(cur_sig(i,:));
            end
        end
    end
end

all_proportions_both_dir = squeeze(mean(all_proportions_both,3,'omitnan'));

%% Calculate proportion increasing
% Output is proportion of changed neurons that are increasing for each
% parameter.
[final_averages] = calc_param_change_direction_24062022(full_all_param_means,all_full_results_both,sig_results_all,sig_summary_cell);
%%

m_means_both = squeeze(mean(all_proportions_both,2,'omitnan'));
m_means_both_dir = squeeze(mean(all_proportions_both_dir,2,'omitnan'));
% titles = ["Left Trials";"Right Trials"];
x_labels = ["Gain/Loss";"No Change";"Location";"Amplitdue";"Width"];

% create increase vs decrease matrix
bar_mat1 = [squeeze(mean(m_means_both_dir(:,1))),squeeze(mean(m_means_both_dir(:,2)));squeeze(mean(m_means_both_dir(:,3))),0]';
bar_mat2 = [squeeze(mean(m_means_both_dir(:,4:end))).*final_averages';squeeze(mean(m_means_both_dir(:,4:end))).*(1-final_averages')];
bar_mat = [bar_mat1,bar_mat2];

bar(bar_mat','stacked')
hold on
for m = 1:num_mice
    scatter([1,2,3,4,5],[m_means_both_dir(m,1)+ m_means_both_dir(m,2),m_means_both_dir(m,3:end)],150,'k','filled')
end
ylim([0,1])
title('Tuning Curve Changes')
ylabel('Proportion Changed')
xticklabels(x_labels)
xtickangle(45)
legend("Increase","Decrease")


%%
figure
% colour_vec = [[0 0.4470 0.7410];[0.8500 0.3250 0.0980];[0.9290 0.6940 0.1250];[0.4940 0.1840 0.5560];[0.4660 0.6740 0.1880];[0.3010 0.7450 0.9330];[0.6350 0.0780 0.1840]];
x_labels = ["Gain";"Loss";"No Change";"Location";"Amplitdue";"Width"];
% Plot totals with no increase/decrease
% bar_mat1 = [squeeze(mean(m_means_both_dir(:,1))),squeeze(mean(m_means_both_dir(:,2)));squeeze(mean(m_means_both_dir(:,3))),0]';
% bar_mat2 = [squeeze(mean(m_means_both_dir(:,4:end))).*final_averages';squeeze(mean(m_means_both_dir(:,4:end))).*(1-final_averages')];
% bar_mat = [bar_mat1,bar_mat2];

bar(squeeze(mean(m_means_both_dir,'omitnan')),'w','LineWidth',2)
hold on
for m = 1:num_mice
    scatter([1,2,3,4,5,6],m_means_both_dir(m,:),150,'k','filled')
end
ylim([0,1])
title('Tuning Curve Changes')
ylabel('Fraction of Neurons')
xticklabels(x_labels)
xtickangle(45)
box off
% axis('square')
% legend("Increase","Decrease")

%% Violin plot version
figure
% colour_vec = [[0 0.4470 0.7410];[0.8500 0.3250 0.0980];[0.9290 0.6940 0.1250];[0.4940 0.1840 0.5560];[0.4660 0.6740 0.1880];[0.3010 0.7450 0.9330];[0.6350 0.0780 0.1840]];
x_labels = ["Gain";"Loss";"No Change";"Location";"Amplitdue";"Width"];
% Plot totals with no increase/decrease
% bar_mat1 = [squeeze(mean(m_means_both_dir(:,1))),squeeze(mean(m_means_both_dir(:,2)));squeeze(mean(m_means_both_dir(:,3))),0]';
% bar_mat2 = [squeeze(mean(m_means_both_dir(:,4:end))).*final_averages';squeeze(mean(m_means_both_dir(:,4:end))).*(1-final_averages')];
% bar_mat = [bar_mat1,bar_mat2];

violin_cell = cell(1,6);
all_proportions_both_dir = permute(all_proportions_both_dir,[3,1,2]);
for i = 1:6
    violin_cell{i} = all_proportions_both_dir(i,:);
end

violin(violin_cell,'facecolor',[0.5,0.5,0.5],'medc',[],'edgecolor',[])
title('Tuning Curve Changes')
ylabel('Fraction of Neurons')
xticklabels(x_labels)
xtickangle(45)
box off
axis('square')
% legend("Increase","Decrease")



%%
% % 
% % figure
% % 
% % subplot(1,2,1)
% % % bar(squeeze(mean(m_means_both_dir)),'w','LineWidth',2)
% % bar(bar_mat','stacked')
% % hold on
% % for m = 1:num_mice
% %     scatter([1,2,3,4,5],[m_means_both_dir(m,1)+ m_means_both_dir(m,2),m_means_both_dir(m,3:end)],150,'k','filled')
% % end
% % ylim([0,1])
% % title('Tuning Curve Changes')
% % ylabel('Proportion Changed')
% % xticklabels(x_labels)
% % xtickangle(30)
% % legend("Increase","Decrease")
% % 

%% Plot 
figure
sig_results = sig_results_all{ex_md(1),ex_md(2)};
% full_results_both = all_full_results_both{ex_md(1),ex_md(2)};
% cur_sig_summary = sig_summary_cell{ex_md(1),ex_md(2)};
cur_sig_check = sig_check_all{ex_md(1),ex_md(2)};
cur_neurons = sum(sig_results==1)>0;

%%
bootstrap_means = bootstrap_means(:,:,:,cur_neurons);
sig_results = sig_results(:,cur_neurons);
% sig_check = sig_check(cur_neurons);
% full_results_both = full_results_both(cur_neurons,:,:);
% cur_sig_summary = cur_sig_summary(:,:,cur_neurons);
cur_sig_check = cur_sig_check(cur_neurons);
zdim = size(sig_results,2);

z_binned_means = zeros(size(bootstrap_means,3),zdim,length(t_types));
  
for i = 1:length(t_types)
    z_binned_means(:,:,i) = mean(squeeze(bootstrap_means(i,:,:,:)),'omitnan');
end

% CIs, upper/lower x ntypes x nbins x zdim
CIs_all = zeros(2,length(t_types),size(bootstrap_means,3),zdim);
for i = 1:size(bootstrap_means,3)
    for j = 1:length(t_types)
        CIs_all(:,j,i,:) = prctile(squeeze(bootstrap_means(j,:,i,:)),CI_vals);
    end
end

% Plot example neurons
titles = ["Gain";"Loss";"No Change";"Location";"Amplitude";"Width"];
num_examples = size(ex_neu,2);
for i = 1:6

    for j = 1:num_examples
        subplot(2,6,(j-1)*6+i)

        cur_n = ex_neu(i,j);
        cur_dir = ex_dir(i,j); % 1 or 2
        for t = 1:2
            % error bars
            % plot(centres,z_binned_means(:,(p-1)*num_neurons+i,t)+squeeze(std(bootstrap_means(t,:,:,(p-1)*num_neurons+i))),'Linewidth',1,'Color',colours_vec(t,:))
            % plot(centres,squeeze(CIs_all(1,cur_dir + 2*(t-1),:,cur_n)),'Linewidth',1,'Color',colours_vec(t,:))
            hold on
            % plot(centres,squeeze(CIs_all(2,cur_dir + 2*(t-1),:,cur_n)),'Linewidth',1,'Color',colours_vec(t,:))
            % plot(centres,z_binned_means(:,(p-1)*num_neurons+i,t)-squeeze(std(bootstrap_means(t,:,:,(p-1)*num_neurons+i))),'Linewidth',1,'Color',colours_vec(t,:))
            h = fill([centres,fliplr(centres)],[squeeze(CIs_all(1,cur_dir + 2*(t-1),:,cur_n))',fliplr(squeeze(CIs_all(2,cur_dir + 2*(t-1),:,cur_n))')],colours_vec(t,:),'LineStyle','none');
            set(h,'facealpha',.3)
            % Means

            plot(centres,z_binned_means(:,cur_n,cur_dir + 2*(t-1)),'Linewidth',3,'Color',colours_vec(t,:))

            % plot baseline
    %         if size(sig_check,1) == length(t_types)
    %             yline(sig_check(t,(p-1)*num_neurons+i),'--','Color',colours_vec(t,:));
    %         end
        end
        xticks([0,200])
        % yticks([])
        % title(cur_n)
        yline(cur_sig_check(cur_n),'--','LineWidth',2);
    end

    %     if size(sig_check,1) == 1
    %         yline(sig_check((p-1)*num_neurons+i),'--','Color','k');
    %     end
        % title(sig_results(1,(p-1)*num_neurons+i)+ " " +sig_results(2,(p-1)*num_neurons+i)+ " " +sig_results(3,(p-1)*num_neurons+i) + " " +sig_results(4,(p-1)*num_neurons+i))
        % title(full_results((p-1)*num_neurons+i,cur_ind,1)+ " " +full_results((p-1)*num_neurons+i,cur_ind,2)+ " " +full_results((p-1)*num_neurons+i,cur_ind,3))
end

for i = 1:6
    subplot(2,6,i)
    title(titles(i))
    % subplot(2,6,6+i)
    % xlabel("Linearised Position (cm)")
end

subplot(2,6,1)
ylabel("\DeltaF/F")

subplot(2,6,7)
ylabel("\DeltaF/F")
xlabel("Linearised Position (cm)")

%% Hierarchical bootstrap

rng(1);
boot_samps = 1000;
num_trials = 4;

all_centres = nan.*ones(size(all_proportions_both_dir,1),1);
all_sems = nan.*ones(size(all_proportions_both_dir,1),1);

stats_data = nan.*ones(size(all_proportions_both_dir,1),num_mice,num_days);
for m = 1:num_mice

    for i = 1:size(all_proportions_both_dir,1)
        cur_data = squeeze(all_proportions_both_dir(i,m,:));

        stats_data(i,m,1:sum(~isnan(cur_data))) = cur_data(~isnan(cur_data));
    end

end

for i = 1:size(all_proportions_both_dir,1)
    bootstats = get_bootstrapped_equalsamples(squeeze(stats_data(i,:,:)),boot_samps,num_trials,'mean');
    %Get mean and SEM of bootstrapped samples:
    all_sems(i) = std(bootstats);
    all_centres(i) = mean(bootstats);
end

h_boots.all_centres = all_centres;
h_boots.all_sems = all_sems;