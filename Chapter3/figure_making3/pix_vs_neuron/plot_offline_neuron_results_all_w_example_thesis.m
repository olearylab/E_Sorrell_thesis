function [stats_tests,h_boots] = plot_offline_neuron_results_all_w_example_thesis(mice,data_days,reps,training_trials)
% 01/12/2022

% Function for making paper ready figures showing comparison between pixel
% and neuron decoding.

% Only neurons and pixels (removed allROIs)

% How much should I separate the different plots?
plot_rows = 2;

res_path = "../shared_files/offline_results/";
res_path_n = "../shared_files/offline_neuron_results_CNN_16112022/";

% set(0,'DefaultAxesFontSize',14)
line_width = 3;
circle_size = 100;
label_size = 50;

colours_vec = [[0 0.4470 0.7410];[0.8500 0.3250 0.0980];[0.9290 0.6940 0.1250];[0.4940 0.1840 0.5560];[0.4660 0.6740 0.1880]];

% example is from mouse 4
ex_mouse = 4;

num_mice = length(mice);
n_days = zeros(num_mice,1);
for m = 1:num_mice
    n_days(m) = length(data_days{m});
end
num_days = max(n_days);

%% Horziontal

figure

% This might not work for other computer screen? set to macbook
% set(gcf,'position',[1,42,1440,763])
set(gcf,'position',[1,42,1440,508])

% Do I plot pitch and yaw raw, or converted into interpretable values?
plot_nums = [1,7,3,4]; % ypos, VA, pitch, yaw
ex_offsets = [1.4804,1.4844];
titles = ["Y Position", "View Angle", "Forward Velocity", "Angular Velocity"];
y_labels = ["Pixels RMSE (m)", "Pixels RMSE (rad)", "Pixels RMSE (m/s)", "Pixels RMSE (rad/s)"];
x_labels = ["Cells RMSE (m)", "Cells RMSE (rad)", "Cells RMSE (m/s)", "Cells RMSE (rad/s)"];
% might need adjusting
ytick_vec = [[0 1 2];[-2,0,2];[0 0.2 0.4];[-1, 0, 1]];

% Example decoding

% results_struct = importdata(res_path + "res_struct_" + example_m_d + "_rep1_" + training_trials + ".mat");
% t = (1:size(results_struct.xtest(results_struct.test_valid,:),1))/results_struct.model_params.fs;

% subtract offset from pitch and yaw, and convert into velocities
circum = 64;
V = 0.32;
alpha = -50/75*circum/V;
beta = 0.05*circum/V/2.5;
% results_struct.xtest(:,3) = alpha*(results_struct.xtest(:,3)-ex_offsets(1));
% results_struct.xprediction(:,3) = alpha*(results_struct.xprediction(:,3)-ex_offsets(1));
% results_struct.xtest(:,4) = -beta*(results_struct.xtest(:,4)-ex_offsets(2));
% results_struct.xprediction(:,4) = -beta*(results_struct.xprediction(:,4)-ex_offsets(2));

% convert from vu into m
vu_conv = 0.0074;
% results_struct.xtest(:,3) = results_struct.xtest(:,3)*vu_conv;
% results_struct.xprediction(:,3) = results_struct.xprediction(:,3)*vu_conv;
% results_struct.xtest(:,1) = results_struct.xtest(:,1)*vu_conv;
% results_struct.xprediction(:,1) = results_struct.xprediction(:,1)*vu_conv;

% Uncomment and edit if including example decoding

ylabel_vec = ["(m)";"(rad)";"(m/s)";"(rad/s)"];
% 
% for i = 1:4
%     ax = subplot(plot_rows,4,i);
%     
%     max1 = max(results_struct.xtest(:,plot_nums(i)));
%     max2 = max(results_struct.xprediction(:,plot_nums(i)));
%     
%     min1 = min(results_struct.xtest(:,plot_nums(i)));
%     min2 = min(results_struct.xprediction(:,plot_nums(i)));
%     
%     baseline = min([min1,min2]);
%     max_val = max([max1,max2]);
%     
%     if i ==2
%         max_val = pi;
%         baseline = -pi;
%     end
%     
%     % Specificly to make slightly nicer for DW113 20210701 example
%     if i ==4
%         max_val = 1.5;
%         baseline = -1.5;
%     end
%     
%     % area(t,(~results_struct.test_valid)*(max_val-baseline)+baseline,baseline,'FaceColor','k','FaceAlpha',0.1,'LineStyle','none')
%     hold on
%     plot(t,results_struct.xtest(results_struct.test_valid,plot_nums(i)),'LineWidth',line_width,'Color',[0.5,0.5,0.6])
%     % plot(t,results_struct.xtest(:,plot_nums(i)),'LineWidth',line_width,'Color',colours_vec(1,:))
%     plot(t,results_struct.xprediction(results_struct.test_valid,plot_nums(i)),'LineWidth',line_width,'Color',colours_vec(ex_mouse,:))
%     title(Titles(i))
%     % ylabel(y_labels(i))
%     xlabel("Time (s)")
%     xlim([0,60])
%     ylim([baseline,max_val]);
%     
%     xticks([0, 60])
%     xticklabels({'0','60'})
%     yticks(ytick_vec(i,:))
%     %yticklabels({'-1','0','1'})
%     ylabel(ylabel_vec(i))
%     
%     box off
%     % axis('square')
% end
% 
% subplot(plot_rows,4,1)
% % text(-17,3,'a','FontSize',label_size,'FontWeight','bold')
% % ylabel("Mouse 4")
% 
% % legend on view angle figure
% subplot(plot_rows,4,2)
% legend('True', 'Decoded')
% 
% subplot(plot_rows,4,4)
% ylim([-1.5,1.5])

% Summary results

% Plot dot and whiskers for each day for each mouse.
% whiskers are full range for now.

colours_vec = [[0 0.4470 0.7410];[0.8500 0.3250 0.0980];[0.9290 0.6940 0.1250];[0.4940 0.1840 0.5560];[0.4660 0.6740 0.1880]];
conv_vec = [vu_conv,1.0,-vu_conv*alpha,-(-beta)];
ytick_vec = [[0 0.25 0.5];[0 0.5 1];[0 0.06 0.12];[0 0.25 0.5]];
ylim_vec = [[0 0.5];[0 1];[0 0.12];[0 0.5]];
%% RMSE
p_rmse = [];
n_rmse = [];
p_rmse_all = [];
n_rmse_all = [];
% r_rmse = [];

% m_means = zeros(length(mice),2);
off_size = 0.3;

m_offset = -off_size*(length(mice)-1)/2:off_size:off_size*(length(mice)-1)/2;

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
            % md_vec_n_all(r,:) = results_struct.all_res(1,plot_nums);
        end
        for i = 1:4
            md_vec(:,i) = md_vec(:,i).*conv_vec(i);
            % md_vec_n_all(:,i) = md_vec_n_all(:,i).*conv_vec(i);
            md_vec_n(:,i) = md_vec_n(:,i).*conv_vec(i);
        end
        m_vec(d,:) = mean(md_vec);
        m_vec_n(d,:) = mean(md_vec_n);
        % m_vec_n_all(d,:) = mean(md_vec_n_all);
        
       
        % Separate sessions
%         p_rmse = [p_rmse;m_vec];
%         n_rmse = [n_rmse;m_vec_n];
%         r_rmse = [r_rmse;m_vec_n_all];
%         for i = 1:4
%             ax = subplot(1,4,i); % edit if including example decoding
%             hold on
%             % Might like to change this for how to plot. Do mean and range
%             % for each mouse. Probably want lines linknig individual mice
%             % as well.
%             scatter((m-1)*3+1,mean(md_vec(:,i)),circle_size,colours_vec(i+1,:),'filled')
%             scatter((m-1)*3+2,mean(md_vec_n_all(:,i)),circle_size,colours_vec(i+1,:),'filled')
%             scatter((m-1)*3+3,mean(md_vec_n(:,i)),circle_size,colours_vec(i+1,:),'filled')
%             
%             ylim([0,1])
%             %errorbar(m,mean(md_vec(:,i)),mean(md_vec(:,i))-min(md_vec(:,i)),max(md_vec(:,i))-mean(md_vec(:,i)),'Color','k')
%             xlim([0,3*length(mice)+1])
%             
%             xticks(1:3*length(mice))
%             xticklabels({'Pixels','ROIs','Cells'}) % Need to edit for more mice
%             yticks([0 0.5 1])
%             yticklabels({'0','0.5','1'})
%             xtickangle(45)
%             title(Titles(i))
%             % xlabel("Mouse")
%         end
    end
    
    p_rmse = [p_rmse;mean(m_vec,1,'omitnan')];
    n_rmse = [n_rmse;mean(m_vec_n,1,'omitnan')];
    % r_rmse = [r_rmse;mean(m_vec_n_all,1,'omitnan')];
    p_rmse_all = [p_rmse_all;m_vec];
    n_rmse_all = [n_rmse_all;m_vec_n];

%% RMSE Plotting removed    
     for i = 1:4
         ax = subplot(plot_rows,4,i); % edit if including example decoding
         hold on
%         % Might like to change this for how to plot. Do mean and range
%         % for each mouse. Probably want lines linknig individual mice
%         % as well.
%         %scatter((m-1)*3+1,mean(md_vec(:,i)),circle_size,colours_vec(i+1,:),'filled')
%         %scatter((m-1)*3+2,mean(md_vec_n_all(:,i)),circle_size,colours_vec(i+1,:),'filled')
%         %scatter((m-1)*3+3,mean(md_vec_n(:,i)),circle_size,colours_vec(i+1,:),'filled')
%         
% %         scatter(1,mean(m_vec(:,i)),circle_size,colours_vec(i+1,:),'filled')
% %         scatter(2,mean(m_vec_n_all(:,i)),circle_size,colours_vec(i+1,:),'filled')
% %         scatter(3,mean(m_vec_n(:,i)),circle_size,colours_vec(i+1,:),'filled')
%      
% % These 2 were the old plots
%      %   scatter(ones(length(data_d),1)+m_offset(m),m_vec(:,i),circle_size,'k','filled')
%      %   scatter(3*ones(length(data_d),1)+m_offset(m),m_vec_n(:,i),circle_size,'k','filled')
%         
%         % errorbar(1,mean(m_vec(:,i)),mean(m_vec(:,i))-min(m_vec(:,i)),max(m_vec(:,i))-mean(m_vec(:,i)),'o','MarkerSize',10,...
%         %            'MarkerEdgeColor',colours_vec(i+1,:),'MarkerFaceColor',colours_vec(i+1,:),'Color','k')
%         % errorbar(2,mean(m_vec_n_all(:,i)),mean(m_vec_n_all(:,i))-min(m_vec_n_all(:,i)),max(m_vec_n_all(:,i))-mean(m_vec_n_all(:,i)),'s','MarkerSize',10,...
%         %            'MarkerEdgeColor',colours_vec(i+1,:),'MarkerFaceColor',colours_vec(i+1,:),'Color','k')
%         % errorbar(3,mean(m_vec_n(:,i)),mean(m_vec_n(:,i))-min(m_vec_n(:,i)),max(m_vec_n(:,i))-mean(m_vec_n(:,i)),'^','MarkerSize',10,...
%         %            'MarkerEdgeColor',colours_vec(i+1,:),'MarkerFaceColor',colours_vec(i+1,:),'Color','k')
%         
%         %% 07032023: scatter cells against pixels
         scatter(m_vec_n(:,i),m_vec(:,i),circle_size,'k','filled')
%         
%         %%
%         
% 
%         % ylim([0,1])
%         %errorbar(m,mean(md_vec(:,i)),mean(md_vec(:,i))-min(md_vec(:,i)),max(md_vec(:,i))-mean(md_vec(:,i)),'Color','k')
%         % plot([1,2,3],[mean(m_vec(:,i)),mean(m_vec_n_all(:,i)),mean(m_vec_n(:,i))],'Color','k')
%         % xlim([0,3*length(mice)+1])
%         % xlim([0,4])
%         % xticks([1,3])
%         % xticks(1:3*length(mice))
%         % xticklabels({'Pixels','Cells'}) 
% 
         yticks(ytick_vec(i,:))
         ylim(ylim_vec(i,:))
         xticks(ytick_vec(i,:))
         xlim(ylim_vec(i,:))
%         % xtickangle(45)
%         % xticks([])
%         % title(Titles(i))
%         % xlabel("Mouse")
         xlabel(x_labels(i))
%         plot(ylim_vec(i,:),ylim_vec(i,:),'--','color','k')
         ylabel(y_labels(i))
     end
end
for i = 1:4
    subplot(plot_rows,4,i)
    plot(ylim_vec(i,:),ylim_vec(i,:),'--','color','k')
    axis('square')
    title(titles(i))
end

% for i = 1:4
%     subplot(plot_rows,4,4+i)
%     % b = bar([1,3],[mean(p_rmse(:,i)),mean(n_rmse(:,i))],'w','LineWidth',2);
%     p1 = plot([0.3,1.7],[mean(p_rmse(:,i)),mean(p_rmse(:,i))],'k','LineWidth',2);
%     uistack(p1,'bottom')
%     p2 = plot([2.3,3.7],[mean(n_rmse(:,i)),mean(n_rmse(:,i))],'k','LineWidth',2);
%     uistack(p2,'bottom')
%     
%     ylabel(ylabel_vec(i))
% end

% subplot(plot_rows,4,4+1)
% ylabel(["Pixels RMSE";"(m)"])


%% R^2
p_R2 = [];
n_R2 = [];

p_R2_all = [];
n_R2_all = [];

p_R2_all_mat = nan.*ones(num_mice,num_days,4);
n_R2_all_mat = nan.*ones(num_mice,num_days,4);
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
            md_vec(r,:) = results_struct.all_res(2,plot_nums);
            results_struct = importdata(res_path_n + "res_struct_" + mice(m) + "_" + data_d(d) + "_rep"+ r + "_" + training_trials + ".mat");
            md_vec_n(r,:) = results_struct.all_res(2,plot_nums);
            % results_struct = importdata(res_path_n + "res_struct_" + mice(m) + "_" + data_d(d) + "_rep"+ r + "_allROI_" + training_trials + ".mat");
            % md_vec_n_all(r,:) = results_struct.all_res(2,plot_nums);
        end
        m_vec(d,:) = mean(md_vec);
        m_vec_n(d,:) = mean(md_vec_n);
        % m_vec_n_all(d,:) = mean(md_vec_n_all);
        
        % Separate essions
%         p_R2 = [p_R2;m_vec];
%         n_R2 = [n_R2;m_vec_n];
%         r_R2 = [r_R2;m_vec_n_all];

%         for i = 1:4
%             ax = subplot(1,4,i); % edit if including example decoding
%             hold on
%             % Might like to change this for how to plot. Do mean and range
%             % for each mouse. Probably want lines linknig individual mice
%             % as well.
%             scatter((m-1)*3+1,mean(md_vec(:,i)),circle_size,colours_vec(i+1,:),'filled')
%             scatter((m-1)*3+2,mean(md_vec_n_all(:,i)),circle_size,colours_vec(i+1,:),'filled')
%             scatter((m-1)*3+3,mean(md_vec_n(:,i)),circle_size,colours_vec(i+1,:),'filled')
%             
%             ylim([0,1])
%             %errorbar(m,mean(md_vec(:,i)),mean(md_vec(:,i))-min(md_vec(:,i)),max(md_vec(:,i))-mean(md_vec(:,i)),'Color','k')
%             xlim([0,3*length(mice)+1])
%             
%             xticks(1:3*length(mice))
%             xticklabels({'Pixels','ROIs','Cells'}) % Need to edit for more mice
%             yticks([0 0.5 1])
%             yticklabels({'0','0.5','1'})
%             xtickangle(45)
%             title(Titles(i))
%             % xlabel("Mouse")
%         end
    end
    % Mice averages
    p_R2 = [p_R2;mean(m_vec,1,'omitnan')];
    n_R2 = [n_R2;mean(m_vec_n,1,'omitnan')];
    
    % all_res
    p_R2_all = [p_R2_all;m_vec];
    n_R2_all = [n_R2_all;m_vec_n];
    % r_R2 = [r_R2;mean(m_vec_n_all,1,'omitnan')];
    
    p_R2_all_mat(m,1:size(m_vec,1),:) = m_vec;
    n_R2_all_mat(m,1:size(m_vec_n,1),:) = m_vec_n;
    
    for i = 1:4
        ax = subplot(plot_rows,4,8+i-4); % edit if including example decoding
        hold on
        % Might like to change this for how to plot. Do mean and range
        % for each mouse. Probably want lines linknig individual mice
        % as well.
        %scatter((m-1)*3+1,mean(md_vec(:,i)),circle_size,colours_vec(i+1,:),'filled')
        %scatter((m-1)*3+2,mean(md_vec_n_all(:,i)),circle_size,colours_vec(i+1,:),'filled')
        %scatter((m-1)*3+3,mean(md_vec_n(:,i)),circle_size,colours_vec(i+1,:),'filled')
        
%         scatter(1,mean(m_vec(:,i)),circle_size,colours_vec(i+1,:),'filled')
%         scatter(2,mean(m_vec_n_all(:,i)),circle_size,colours_vec(i+1,:),'filled')
%         scatter(3,mean(m_vec_n(:,i)),circle_size,colours_vec(i+1,:),'filled')
        
      %  scatter(ones(length(data_d),1)+m_offset(m),m_vec(:,i),circle_size,'k','filled')
      %  scatter(3*ones(length(data_d),1)+m_offset(m),m_vec_n(:,i),circle_size,'k','filled')
        
        % errorbar(1,mean(m_vec(:,i)),mean(m_vec(:,i))-min(m_vec(:,i)),max(m_vec(:,i))-mean(m_vec(:,i)),'o','MarkerSize',10,...
        %             'MarkerEdgeColor',colours_vec(i+1,:),'MarkerFaceColor',colours_vec(i+1,:),'Color','k')
        % errorbar(2,mean(m_vec_n_all(:,i)),mean(m_vec_n_all(:,i))-min(m_vec_n_all(:,i)),max(m_vec_n_all(:,i))-mean(m_vec_n_all(:,i)),'s','MarkerSize',10,...
        %             'MarkerEdgeColor',colours_vec(i+1,:),'MarkerFaceColor',colours_vec(i+1,:),'Color','k')
        % errorbar(3,mean(m_vec_n(:,i)),mean(m_vec_n(:,i))-min(m_vec_n(:,i)),max(m_vec_n(:,i))-mean(m_vec_n(:,i)),'^','MarkerSize',10,...
        %             'MarkerEdgeColor',colours_vec(i+1,:),'MarkerFaceColor',colours_vec(i+1,:),'Color','k')
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
        yticks([0 0.5 1])
        yticklabels({'0','0.5','1'})
        xticks([0 0.5 1])
        xticklabels({'0','0.5','1'})
        % xtickangle(45)
        % title(Titles(i))
        % xlabel("Mouse")
        
        xlabel("Cells R^{2}")
    end
end
for i = 1:4
    subplot(plot_rows,4,8+i-4)
    plot([0,1],[0,1],'--','color','k')
    axis('square')
    ylabel("Pixels R^{2}")
end

% for i = 1:4
%     subplot(plot_rows,4,8+i)
%     % b = bar([1,3],[mean(p_R2(:,i)),mean(n_R2(:,i))],'w','LineWidth',2);
%     p1 = plot([0.3,1.7],[mean(p_R2(:,i)),mean(p_R2(:,i))],'k','LineWidth',2);
%     uistack(p1,'bottom')
%     p2 = plot([2.3,3.7],[mean(n_R2(:,i)),mean(n_R2(:,i))],'k','LineWidth',2);
%     uistack(p2,'bottom')
% end

% subplot(plot_rows,4,8+1-4)
% ylabel("Pixels R^{2}")

% text(-1.7,1.1,'b','FontSize',label_size,'FontWeight','bold')

%% Statistical tests
% Perform two sided-Wilcoxon signed rank test, and two sided paired sample
% t-test, using each mouse as a data-point (14 sessions for 5 mice)
pn_wilc = zeros(2,4);
% pr_wilc = zeros(2,4);

pn_t = zeros(2,4);
% pr_t = zeros(2,4);

pn_wilc_all = zeros(2,4);

p_medians_all = zeros(2,4);
n_medians_all = zeros(2,4);

for i = 1:4
    pn_wilc(1,i) = signrank(p_rmse(:,i),n_rmse(:,i));
    pn_wilc(2,i) = signrank(p_R2(:,i),n_R2(:,i));
    
    pn_wilc_all(1,i) = signrank(p_rmse_all(:,i),n_rmse_all(:,i));
    pn_wilc_all(2,i) = signrank(p_R2_all(:,i),n_R2_all(:,i));
    
    p_medians_all(1,i) = median(p_rmse(:,i));
    p_medians_all(2,i) = median(p_R2(:,i));
    n_medians_all(1,i) = median(n_rmse(:,i));
    n_medians_all(2,i) = median(n_R2(:,i));
    
    % pr_wilc(1,i) = signrank(p_rmse(:,i),r_rmse(:,i));
    % pr_wilc(2,i) = signrank(p_R2(:,i),r_R2(:,i));
    
    [~,pn_t(1,i)] = ttest(p_rmse(:,i),n_rmse(:,i));
    [~,pn_t(2,i)] = ttest(p_R2(:,i),n_R2(:,i));
    
    % [~,pr_t(1,i)] = ttest(p_rmse(:,i),r_rmse(:,i));
    % [~,pr_t(2,i)] = ttest(p_R2(:,i),r_R2(:,i));
end
stats_tests.pn_wilc = pn_wilc;
stats_tests.pn_wilc_all = pn_wilc_all;
stats_tests.p_medians_all = p_medians_all;
stats_tests.n_medians_all = n_medians_all;
% stats_tests.pr_wilc = pr_wilc;
stats_tests.pn_t = pn_t;
% stats_tests.pr_t = pr_t;

%% Hierarchical Bootstrap

% now for neuron pixel comparison
all_p_boot = nan.*ones(4,1);
all_centres = nan.*ones(4,2);
all_sems = nan.*ones(4,2);

for i = 1:4
    [all_p_boot(i),all_centres(i,:),all_sems(i,:)] = run_H_boot_ets(squeeze(p_R2_all_mat(:,:,i)), squeeze(n_R2_all_mat(:,:,i)),false);
end

h_boots.all_p_boot = all_p_boot;
h_boots.all_centres = all_centres;
h_boots.all_sems = all_sems;
