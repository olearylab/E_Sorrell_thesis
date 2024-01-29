function [h_boots] = plot_offline_decoder_comparisons(mice,data_days,example_m_d,reps,training_trials)
% 06/06/2023

% Function for making paper ready figures showing comparison between
% different decoding algorithms

% Plot example decoding for each of the algorithms, then comparison of
% accuracy (RMSE)
% plot_rows = 5;

res_path = "../shared_files/offline_results/";
res_path_k = "offline_results_kalman/";
% res_path_b = "offline_results_bayes/";
res_path_b = "offline_results_bayes_alt_prior/";
res_path_p = "offline_results_particle/";
all_res_paths = [res_path;res_path_k;res_path_b;res_path_p];

% set(0,'DefaultAxesFontSize',14)
line_width = 3;
circle_size = 100;
label_size = 50;

colours_vec = [[0 0.4470 0.7410];[0.8500 0.3250 0.0980];[0.9290 0.6940 0.1250];[0.4940 0.1840 0.5560];[0.4660 0.6740 0.1880]];

% example is from mouse 4
ex_mouse = 4;

%% Horziontal

% figure

% This might not work for other computer screen? set to macbook
% set(gcf,'position',[1,42,1440,763])
% set(gcf,'position',[1,42,1440,508])

% Do I plot pitch and yaw raw, or converted into interpretable values?
plot_nums = [1,7,3,4;1,2,3,4;1,2,3,4;1,2,3,4]; % ypos, VA, pitch, yaw
ex_offsets = [1.4804,1.4844];
Titles = ["Y Position", "View Angle", "Forward Velocity", "Angular Velocity"];
% y_labels = ["Y Position (m)", "View Angle (rad)", "Forward Velocity (m/s)", "Angular Velocity (rad/s)"];
% might need adjusting
ytick_vec = [[0 1 2];[-2,0,2];[0 0.2 0.4];[-1, 0, 1]];

% Example decoding
for ii = 1:4
    figure
    set(gcf,'position',[1,42,1440,254])
    results_struct = importdata(all_res_paths(ii) + "res_struct_" + example_m_d + "_rep1_" + training_trials + ".mat");
    t = (1:size(results_struct.xtest(results_struct.test_valid,:),1))/results_struct.model_params.fs;

    % subtract offset from pitch and yaw, and convert into velocities
    circum = 64;
    V = 0.32;
    alpha = -50/75*circum/V;
    beta = 0.05*circum/V/2.5;
    results_struct.xtest(:,3) = alpha*(results_struct.xtest(:,3)-ex_offsets(1));
    results_struct.xprediction(:,3) = alpha*(results_struct.xprediction(:,3)-ex_offsets(1));
    results_struct.xtest(:,4) = -beta*(results_struct.xtest(:,4)-ex_offsets(2));
    results_struct.xprediction(:,4) = -beta*(results_struct.xprediction(:,4)-ex_offsets(2));

    % convert from vu into m
    vu_conv = 0.0074;
    results_struct.xtest(:,3) = results_struct.xtest(:,3)*vu_conv;
    results_struct.xprediction(:,3) = results_struct.xprediction(:,3)*vu_conv;
    results_struct.xtest(:,1) = results_struct.xtest(:,1)*vu_conv;
    results_struct.xprediction(:,1) = results_struct.xprediction(:,1)*vu_conv;

    % Uncomment and edit if including example decoding

    ylabel_vec = ["(m)";"(rad)";"(m/s)";"(rad/s)"];

    for i = 1:4
        ax = subplot(1,4,i);

        max1 = max(results_struct.xtest(:,plot_nums(ii,i)));
        max2 = max(results_struct.xprediction(:,plot_nums(ii,i)));

        min1 = min(results_struct.xtest(:,plot_nums(ii,i)));
        min2 = min(results_struct.xprediction(:,plot_nums(ii,i)));

        baseline = min([min1,min2]);
        max_val = max([max1,max2]);

        if i ==2
            max_val = pi;
            baseline = -pi;
        end

        % Specificly to make slightly nicer for DW113 20210701 example
        if i ==4
            max_val = 1.5;
            baseline = -1.5;
        end

        % area(t,(~results_struct.test_valid)*(max_val-baseline)+baseline,baseline,'FaceColor','k','FaceAlpha',0.1,'LineStyle','none')
        hold on
        plot(t,results_struct.xtest(results_struct.test_valid,plot_nums(ii,i)),'LineWidth',line_width,'Color',[0.5,0.5,0.6])
        % plot(t,results_struct.xtest(:,plot_nums(i)),'LineWidth',line_width,'Color',colours_vec(1,:))
        plot(t,results_struct.xprediction(results_struct.test_valid,plot_nums(ii,i)),'LineWidth',line_width,'Color',colours_vec(ex_mouse,:))
        title(Titles(i))
        % ylabel(y_labels(i))
        xlabel("Time (s)")
        xlim([0,60])
        ylim([baseline,max_val]);

        xticks([0, 60])
        xticklabels({'0','60'})
        yticks(ytick_vec(i,:))
        %yticklabels({'-1','0','1'})
        ylabel(ylabel_vec(i))

        box off
        % axis('square')
    end

    % subplot(plot_rows,4,1)
    % text(-17,3,'a','FontSize',label_size,'FontWeight','bold')
    % ylabel("Mouse 4")

    % legend on view angle figure
    if ii == 1
        subplot(1,4,2)
        legend('True', 'Decoded')
    end

    subplot(1,4,4)
    ylim([-1.5,1.5])
    
end

% Summary results

% Violin plot comparisons

colours_vec = [[0 0.4470 0.7410];[0.8500 0.3250 0.0980];[0.9290 0.6940 0.1250];[0.4940 0.1840 0.5560];[0.4660 0.6740 0.1880]];
conv_vec = [vu_conv,1.0,-vu_conv*alpha,-(-beta)];
ytick_vec = [[0 0.25 0.5];[0 0.5 1];[0 0.06 0.12];[0 0.25 0.5]];
ylim_vec = [[0 0.5];[0 1];[0 0.12];[0 0.5]];
%% RMSE

% m_means = zeros(length(mice),2);
off_size = 0.3;

m_offset = -off_size*(length(mice)-1)/2:off_size:off_size*(length(mice)-1)/2;

v_cell = cell(1,4);
cur_cell = cell(4,4);

num_mice = length(mice);
n_days = zeros(num_mice,1);
for m = 1:num_mice
    n_days(m) = length(data_days{m});
end
num_days = max(n_days);

data_cell = cell(1,4);

for ii = 1:4
    p_rmse = [];
    p_rmse_all = [];

    p_rmse_all_mat = nan.*ones(num_mice,num_days,4);
    
    for m = 1:length(mice)
        data_d = data_days{m};
        m_vec = zeros(length(data_d),4);
        for d = 1:length(data_d)
            md_vec = zeros(reps,4);
            for r = 1:reps
                results_struct = importdata(all_res_paths(ii) + "res_struct_" + mice(m) + "_" + data_d(d) + "_rep"+ r + "_" + training_trials + ".mat");
                md_vec(r,:) = results_struct.all_res(1,plot_nums(ii,:));
            end
            for i = 1:4
                md_vec(:,i) = md_vec(:,i).*conv_vec(i);
            end
            m_vec(d,:) = mean(md_vec);            
        end

        p_rmse = [p_rmse;mean(m_vec,1,'omitnan')];
        p_rmse_all = [p_rmse_all;m_vec];
        
        p_rmse_all_mat(m,1:size(m_vec,1),:) = m_vec;

    end
    for i = 1:4
        cur_cell{i,ii} = p_rmse_all(:,i);
    end
    data_cell{ii} = p_rmse_all_mat;
end

for i = 1:4
    temp_cell = cell(1,4);
    for ii = 1:4
        temp_cell{ii} = cur_cell{i,ii};
    end
    v_cell{i} = temp_cell;
end

% Plot violin
y_labs = ["RMSE (m)";"RMSE (rad)";"RMSE (m/s)";"RMSE (rad/s)"];
figure
set(gcf,'position',[1,42,1440,254])
for i = 1:4
    subplot(1,4,i)
    violin(v_cell{i},'facecolor',[0.5,0.5,0.5],'medc',[],'edgecolor',[])
    box off
    ylabel(y_labs(i))
    xticks([1,2,3,4])
    xticklabels(["Linear";"Kalman";"MAP";"Particle"])
    xtickangle(30)
    legend off
end

%% Hierarchical Bootstrap

% now for neuron pixel comparison
all_p_boot = nan.*ones(4,4);
all_centres = nan.*ones(4,4);
all_sems = nan.*ones(4,4);

for i = 1:4
    for j = 1:4
        [all_p_boot(i,j),h_centres,h_sems] = run_H_boot_ets(squeeze(data_cell{i}(:,:,i)), squeeze(data_cell{j}(:,:,i)),false);
        all_centres(i,j) = h_centres(2);
        all_sems(i,j) = h_sems(2);
    end
end

h_boots.all_p_boot = all_p_boot;
h_boots.all_centres = all_centres;
h_boots.all_sems = all_sems;
