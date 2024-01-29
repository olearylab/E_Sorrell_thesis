function [h_boots] = plot_preprocessing_quant(mice,data_days,reps,training_trials)
% 29/05/2023

res_path = "offline_results/";
res_path_n = "offline_results_raw/";

line_width = 3;
circle_size = 100;

vu_conv = 0.0074;
plot_nums = [1,7,3,4]; % ypos, VA, pitch, yaw
circum = 64;
V = 0.32;
alpha = -50/75*circum/V;
beta = 0.05*circum/V/2.5;
titles = ["Y Position (m)", "View Angle (rad)", "Forward Velocity (m/s)", "Angular Velocity (rad/s)"];
colours_vec = [[0 0.4470 0.7410];[0.8500 0.3250 0.0980];[0.9290 0.6940 0.1250];[0.4940 0.1840 0.5560];[0.4660 0.6740 0.1880]];
conv_vec = [vu_conv,1.0,-vu_conv*alpha,-(-beta)];
ytick_vec = [[0 0.5 1];[0 0.5 1];[0 0.06 0.12];[0 0.25 0.5]];
ylim_vec = [[0 1];[0 1];[0 0.12];[0 0.5]];

figure 
%% R^2
p_rmse = [];
n_rmse = [];
p_rmse_all = [];
n_rmse_all = [];

num_mice = length(mice);
n_days = zeros(num_mice,1);
for m = 1:num_mice
    n_days(m) = length(data_days{m});
end
num_days = max(n_days);

p_rmse_all_mat = nan.*ones(num_mice,num_days,4);
n_rmse_all_mat = nan.*ones(num_mice,num_days,4);
% r_R2 = [];
% Need to edit if going to actually include multiple mice
for m = 1:length(mice)
    data_d = data_days{m};
    m_vec = zeros(length(data_d),4);
    m_vec_n = zeros(length(data_d),4);
    % m_vec_n_all = zeros(length(data_d),4);
    for d = 1:length(data_d)
        md_vec = zeros(reps,4);
        md_vec_n = zeros(reps,4);
        % md_vec_n_all = zeros(reps,4);
        for r = 1:reps
            results_struct = importdata(res_path + "res_struct_" + mice(m) + "_" + data_d(d) + "_rep"+ r + "_" + training_trials + ".mat");
            md_vec(r,:) = results_struct.all_res(1,plot_nums);
            results_struct = importdata(res_path_n + "res_struct_" + mice(m) + "_" + data_d(d) + "_rep"+ r + "_" + training_trials + ".mat");
            md_vec_n(r,:) = results_struct.all_res(1,plot_nums);
            % results_struct = importdata(res_path_n + "res_struct_" + mice(m) + "_" + data_d(d) + "_rep"+ r + "_allROI_" + training_trials + ".mat");
            % md_vec_n_all(r,:) = results_struct.all_res(2,plot_nums);
        end
        for i = 1:4
            md_vec(:,i) = md_vec(:,i).*conv_vec(i);
            % md_vec_n_all(:,i) = md_vec_n_all(:,i).*conv_vec(i);
            md_vec_n(:,i) = md_vec_n(:,i).*conv_vec(i);
        end
        m_vec(d,:) = mean(md_vec);
        m_vec_n(d,:) = mean(md_vec_n);
        % m_vec_n_all(d,:) = mean(md_vec_n_all);
    end
    % Mice averages
    p_rmse = [p_rmse;mean(m_vec,1,'omitnan')];
    n_rmse = [n_rmse;mean(m_vec_n,1,'omitnan')];
    % r_rmse = [r_rmse;mean(m_vec_n_all,1,'omitnan')];
    p_rmse_all = [p_rmse_all;m_vec];
    n_rmse_all = [n_rmse_all;m_vec_n];
    
    p_rmse_all_mat(m,1:size(m_vec,1),:) = m_vec;
    n_rmse_all_mat(m,1:size(m_vec_n,1),:) = m_vec_n;
    for i = 1:4
        ax = subplot(1,4,i); % edit if including example decoding
        hold on
        
        % Might like to do violin or something instead of scatter?
%% 07032023: scatter cells against pixels
        scatter(m_vec_n(:,i),m_vec(:,i),circle_size,'k','filled')
        
        %%
        % ylim([0,1])
        %errorbar(m,mean(md_vec(:,i)),mean(md_vec(:,i))-min(md_vec(:,i)),max(md_vec(:,i))-mean(md_vec(:,i)),'Color','k')
        % plot([1,2,3],[mean(m_vec(:,i)),mean(m_vec_n_all(:,i)),mean(m_vec_n(:,i))],'Color','k')
        % xlim([0,3*length(mice)+1])
        % xlim([0,4])
        % xticks([1,3])
        % xticks(1:3*length(mice))
        % xticklabels({'Pixels','Cells'}) % Need to edit for more mice
        yticks(ytick_vec(i,:))
        ylim(ylim_vec(i,:))
        xticks(ytick_vec(i,:))
        xlim(ylim_vec(i,:))
        % xtickangle(45)
        % title(Titles(i))
        % xlabel("Mouse")
        axis('equal')
        xlabel("Raw RMSE")
    end
end
% ylim_vec = [[0 1];[0 1];[0 0.12];[0 0.5]];
for i = 1:4
    subplot(1,4,i)
    title(titles(i))
    plot(ylim_vec(i,:),ylim_vec(i,:),'--','color','k')
    axis('square')
end

% for i = 1:4
%     subplot(plot_rows,4,8+i)
%     % b = bar([1,3],[mean(p_R2(:,i)),mean(n_R2(:,i))],'w','LineWidth',2);
%     p1 = plot([0.3,1.7],[mean(p_R2(:,i)),mean(p_R2(:,i))],'k','LineWidth',2);
%     uistack(p1,'bottom')
%     p2 = plot([2.3,3.7],[mean(n_R2(:,i)),mean(n_R2(:,i))],'k','LineWidth',2);
%     uistack(p2,'bottom')
% end

subplot(1,4,1)
ylabel("Preprocessed RMSE")

%% Hierarchical Bootstrap

% now for neuron pixel comparison
all_p_boot = nan.*ones(4,1);
all_centres = nan.*ones(4,2);
all_sems = nan.*ones(4,2);

for i = 1:4
    [all_p_boot(i),all_centres(i,:),all_sems(i,:)] = run_H_boot_ets(squeeze(p_rmse_all_mat(:,:,i)), squeeze(n_rmse_all_mat(:,:,i)),false);
end

h_boots.all_p_boot = all_p_boot;
h_boots.all_centres = all_centres;
h_boots.all_sems = all_sems;

