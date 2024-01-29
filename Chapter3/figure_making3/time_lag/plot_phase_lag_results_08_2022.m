function [] = plot_phase_lag_results_08_2022(all_results,lag_vec)
% 19/08/2022

% Plot results of decoding at different phase shifts.
% Plot for many mice and days
mice_vec = [1,1,1,2,2,2,3,3,3,4,4,4,4,5];
day_vec = [1,2,3,1,2,3,1,2,3,1,2,3,4,1];
sess_mat = NaN(length(lag_vec),4,5,4);

num_sess = length(all_results);
num_lags = length(lag_vec);

mid_ind = ceil(length(lag_vec)/2);

figure

% colour_vec = [[0 0.4470 0.7410];[0.8500 0.3250 0.0980];[0.9290 0.6940 0.1250];[0.4940 0.1840 0.5560];[0.4660 0.6740 0.1880];[0.3010 0.7450 0.9330];[0.6350 0.0780 0.1840]];
all_sess_res = nan.*ones(num_sess,num_lags,4);
for m = 1:num_sess
    cur_res = all_results{m};

    % percent change in RMSE for ypos, va, pitch and yaw
    % Do I need to convert into correct units before calculating? Am not for
    % now.
    plot_titles = ["Y position","View angle","Forward velocity","Angular velocity"];
    plot_inds = [1,7,3,4];
   
    cur_plot = zeros(num_lags,10);
    for l = 1:num_lags
        cur_all_res = cur_res{l}.all_res;
        cur_plot(l,:) = cur_all_res(1,:);
    end
    for j = 1:4
        subplot(1,4,j)
        plot(lag_vec./30,100.*((cur_plot(:,plot_inds(j))-cur_plot(mid_ind,plot_inds(j)))./cur_plot(mid_ind,plot_inds(j))),'LineWidth',2,'Color',[0.5,0.5,0.5])
        hold on
        all_sess_res(m,:,j) = 100.*((cur_plot(:,plot_inds(j))-cur_plot(mid_ind,plot_inds(j)))./cur_plot(mid_ind,plot_inds(j)));
        sess_mat(:,j,mice_vec(m),day_vec(m)) = 100.*((cur_plot(:,plot_inds(j))-cur_plot(mid_ind,plot_inds(j)))./cur_plot(mid_ind,plot_inds(j)));
    end
    
    % sess_mat(:,:,mice_vec(m),day_vec(m)) = cur_plot(:,plot_inds);
    
end
for j =1:4
    subplot(1,4,j)
    plot(lag_vec./30,squeeze(mean(all_sess_res(:,:,j),'omitnan')),'LineWidth',3,'Color','k')
    xline(0,'--','LineWidth',2);
    yline(0,'--','LineWidth',2);
    box off
    ylim([-30,40])
    xlabel("Time Lag (s)")
    xticks([-1,0,1])
    % yticks([])
    title(plot_titles(j))
    ylabel("% Change in RMSE")
    axis('square')
end

% subplot(1,4,1)
% ylabel("% Change in RMSE")
% yticks([-20,-10,0,10,20])

%% H boot plot
figure
all_centres = NaN(length(lag_vec),4);
all_sems = NaN(length(lag_vec),4);
all_p_boot = NaN(length(lag_vec),2);
for b = 1:length(lag_vec)
    [all_p_boot(b,1),all_centres(b,1:2),all_sems(b,1:2)] = run_H_boot_ets(squeeze(sess_mat(b,1,:,:)), squeeze(sess_mat(b,2,:,:)),false);
    [all_p_boot(b,2),all_centres(b,3:4),all_sems(b,3:4)] = run_H_boot_ets(squeeze(sess_mat(b,3,:,:)), squeeze(sess_mat(b,4,:,:)),false);
end

h_boots.all_p_boot = all_p_boot;
h_boots.all_centres = all_centres;
h_boots.all_sems = all_sems;

figure

% boostrap stds instead
boot_stds = all_sems;
lims_all = zeros(2,length(lag_vec),4);
lims_all(1,:,:) = all_centres - boot_stds;
lims_all(2,:,:) = all_centres + boot_stds;

for i = 1:4
    subplot(1,4,i)
    % h = fill([centres,fliplr(centres)],[CIs_all(1,:),fliplr(CIs_all(2,:))],'k','EdgeColor','none');
    h = fill([lag_vec./30,fliplr(lag_vec./30)],[squeeze(lims_all(1,:,i)),fliplr(squeeze(lims_all(2,:,i)))],'k','EdgeColor','none');
    set(h,'facealpha',.3)
    hold on

    plot(lag_vec./30,all_centres(:,i),'LineWidth',2,'Color','k')

    % plot(lag_vec./30,squeeze(mean(all_sess_res(:,:,j),'omitnan')),'LineWidth',3,'Color','k')
    xline(0,'--','LineWidth',2);
    yline(0,'--','LineWidth',2);
    box off
    ylim([-20,20])
    xlabel("Time lag (s)")
    xticks([-1,0,1])
    % yticks([])
    title(plot_titles(i))
    ylabel("% Change in RMSE")
    axis('square')

end