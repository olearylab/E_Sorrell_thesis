function [] = plot_pva_trajectory_comparisons_0612022(virmen_data,tbt_details,nbins,linearise_x,boot_samps,sub_sample,ex_offset)
% 0612022

% function for plotting comparison of pitch and view angle between ball and
% bmi trials (i.e. decoder output on bmi trials, to control signal on ball
% trials.

% For adding into figure 3

% Potentially also add a summary figure across all mice and sessions
% (compare across to within type variation.

% Might be hard to plot left and right on same figure, so might need 4
% figures

x_vec = [13,7,16,7];
% Should I use m/s or cm/s
pos_scale = 0.74;
ex_mouse = 4;

circum = 64;
V = 0.32;
alpha = -50/75*circum/V;
beta = 0.05*circum/V/2.5;

virmen_data(13,:) = (alpha*(virmen_data(13,:)-ex_offset))*pos_scale;
virmen_data(16,:) = (alpha*(virmen_data(16,:)-ex_offset))*pos_scale;

virmen_data(7,:) = wrapToPi(virmen_data(7,:));

all_colours = [[0 0.4470 0.7410];[0.8500 0.3250 0.0980];[0.9290 0.6940 0.1250];[0.4940 0.1840 0.5560];[0.4660 0.6740 0.1880]];

%%

% Do I want individual trials plotted? Or mean and std, or CIs?
% Do I need to subsample? Bootstrap?

%% Individual trials
% num trials x nbins x numx
[x_binned,centres] = bin_kin_data(virmen_data,x_vec,linearise_x,nbins);
centres = centres*pos_scale;

types_vec = [1,4,7,10];
line_vec = ["-";"--";"-";"--"];
% color_vec = [[0 0 0];[0 0 0];all_colours(ex_mouse,:);all_colours(ex_mouse,:)];
color_vec = [[0.4940 0.1840 0.5560];[0.4940 0.1840 0.5560];[0.4660 0.6740 0.1880];[0.4660 0.6740 0.1880]];

%% Mean and STD
% 
% figure
% for i = 1:2
%     subplot(1,2,i)
%     for j = 1:4
%         cur_mean = squeeze(mean(x_binned(tbt_details(3,:)==types_vec(j),:,i)));
%         cur_std = squeeze(std(x_binned(tbt_details(3,:)==types_vec(j),:,i)));
%         
%         plot(centres,cur_mean-cur_std,line_vec(j),'Color',color_vec(j,:),'LineWidth',1)
%         hold on
%         plot(centres,cur_mean+cur_std,line_vec(j),'Color',color_vec(j,:),'LineWidth',1)
%         h = fill([centres,fliplr(centres)],[(cur_mean-cur_std),fliplr((cur_mean+cur_std))],color_vec(j,:));
%         set(h,'facealpha',.1)
%         
%         plot(centres,cur_mean,line_vec(j),'Color',color_vec(j,:),'LineWidth',2)
%         
%         
%     end
% end


%% Bootstrapped with CIs

% ex_colour = [0 0.4470 0.7410];

% color_vec = [[0 0 0];all_colours(ex_mouse,:)];
% color_vec = [[0 0 0];ex_colour];
color_vec = [[0.4940 0.1840 0.5560];[0.4660 0.6740 0.1880]];

[bootstrap_means] = calc_bootsrapped_means_subsample(x_binned,tbt_details,boot_samps,types_vec,sub_sample);
x_binned_means = zeros(size(bootstrap_means,3),length(x_vec),length(types_vec));
  
for i = 1:length(types_vec)
    x_binned_means(:,:,i) = mean(squeeze(bootstrap_means(i,:,:,:)),'omitnan');
end

% CI_vals = [16,84];
CI_vals = [2.5,97.5];

% CIs, upper/lower x ntypes x nbins x zdim
CIs_all = zeros(2,length(types_vec),size(bootstrap_means,3),length(x_vec));
for i = 1:size(bootstrap_means,3)
    for j = 1:length(types_vec)
        CIs_all(:,j,i,:) = prctile(squeeze(bootstrap_means(j,:,i,:)),CI_vals);
    end
end

figure
for i = 1:2
    for j = 1:2        
        subplot(2,2,(i-1)*2+j)
%         plot(centres,squeeze(CIs_all(1,j,:,i)),'Color',color_vec(1,:),'LineWidth',1)
         hold on
%         plot(centres,squeeze(CIs_all(2,j,:,i)),'Color',color_vec(1,:),'LineWidth',1)
        h = fill([centres,fliplr(centres)],[squeeze(CIs_all(1,j,:,i))',fliplr(squeeze(CIs_all(2,j,:,i))')],color_vec(1,:),'EdgeColor','none');
        set(h,'facealpha',.3)
        
        plot(centres,squeeze(x_binned_means(:,i,j)),'Color',color_vec(1,:),'LineWidth',2)
        
%         plot(centres,squeeze(CIs_all(1,j+2,:,i)),'Color',color_vec(2,:),'LineWidth',1)
         hold on
%         plot(centres,squeeze(CIs_all(2,j+2,:,i)),'Color',color_vec(2,:),'LineWidth',1)
        h = fill([centres,fliplr(centres)],[squeeze(CIs_all(1,j+2,:,i+2))',fliplr(squeeze(CIs_all(2,j+2,:,i+2))')],color_vec(2,:),'EdgeColor','none');
        set(h,'facealpha',.3)
        
        plot(centres,squeeze(x_binned_means(:,i+2,j+2)),'Color',color_vec(2,:),'LineWidth',2)
        
        if i == 2
            ylim([-2,2]);
            yline(0,'--','LineWidth',2);
        else
            % insert ylim for pitch
            ylim([0,40]);
        end
        box off
    end
end
subplot(2,2,1)
ylabel(["Forward velocity"; "(cm/s)"]);
title("Left Trials")
% axis('square')

subplot(2,2,2)
title("Right trials")
% axis('square')

subplot(2,2,3)
ylabel(["View angle"; "(rad)"]);
xlabel("Linearized position (cm)")
% axis('square')

subplot(2,2,4)
xlabel("Linearized position (cm)")
% axis('square')
