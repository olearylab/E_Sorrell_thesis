function [] = plot_example_trajectories_06122022(virmen_cell, tbt_cell, ex_md, offsets, linearise_x, nbins)
% 06/12/2022

% Plot example trajectories for pitch and yaw to compare between bmi and
% control trials.

% Maybe somehow combine with LSTM figures.
x_vec = [13,15];
%% 
% 2 x 2
% Best to collapse into 1 x 2 if possible
% Do I want actual time? Or binned?
types_vec = [1,4,7,10];

virmen_data = virmen_cell{ex_md(1),ex_md(2)};
tbt_details = tbt_cell{ex_md(1),ex_md(2)};

ITI = virmen_data(8,:);
cleaned_valid = clean_valid_data(ITI);
virmen_data = virmen_data(:,cleaned_valid);

% convert into velocity
circum = 64;
V = 0.32;
alpha = -50/75*circum/V;
beta = 0.05*circum/V/2.5;
virmen_data(13,:) = alpha*(virmen_data(13,:)-offsets(1));
virmen_data(15,:) = -beta*(virmen_data(15,:)-offsets(2));

% Bin data: trials x nbins x xdim
[x_binned,centres] = bin_kin_data(virmen_data,x_vec,linearise_x,nbins);

%% Plot averages and error bars - alternative
% Plot either mean or median line with shaded in standard deviation, or
% confidence interval, or something similar

% Pick better colour scheme?
plot_types = [1,7;4,10];
% Do I want lighter indivudal trials? Or lighter stds or something similar?
% Do I need to have separate colours in this case?
% plot_colour = [0 0.4470 0.7410];
plot_colours = [[0.4940 0.1840 0.5560];[0.4660 0.6740 0.1880]];
figure
pos_scale = .74;
centres = centres*pos_scale;

% with std
% for i = 1:2
%     if sum(tbt_details(3,:)==plot_types(i,1)) == 1
%         cur_mean_n = squeeze(x_binned(tbt_details(3,:)==plot_types(i,1),:,:));
%         cur_mean_b = squeeze(x_binned(tbt_details(3,:)==plot_types(i,2),:,:));
%         % maybe change from std
%         cur_std_n = [0,0];
%         cur_std_b = [0,0];
%     else
%         cur_mean_n = squeeze(mean(x_binned(tbt_details(3,:)==plot_types(i,1),:,:),'omitnan'));
%         cur_mean_b = squeeze(mean(x_binned(tbt_details(3,:)==plot_types(i,2),:,:),'omitnan'));
%         cur_std_n = squeeze(std(x_binned(tbt_details(3,:)==plot_types(i,1),:,:),'omitnan'))./sqrt(length(find(tbt_details(3,:)==plot_types(i,1))));
%         cur_std_b = squeeze(std(x_binned(tbt_details(3,:)==plot_types(i,2),:,:),'omitnan'))./sqrt(length(find(tbt_details(3,:)==plot_types(i,2))));
%     end
%     for j = 1:2
%         % subplot(1,4,2*(j-1)+i)
%         subplot(2,2,2*(j-1)+i)
%         yline(0,'LineWidth',2);
%         hold on
%         plot(centres,cur_mean_n(:,j)-cur_std_n(:,j),'Color',plot_colours(1,:),'LineWidth',1)
%         plot(centres,cur_mean_n(:,j)+cur_std_n(:,j),'Color',plot_colours(1,:),'LineWidth',1)
%         % xx = 1:100;
%         h = fill([centres,fliplr(centres)],[(cur_mean_n(:,j)-cur_std_n(:,j))',fliplr((cur_mean_n(:,j)+cur_std_n(:,j))')],plot_colours(1,:));
%         set(h,'facealpha',.1)
%         plot(centres,cur_mean_n(:,j),'LineWidth',2,'Color',plot_colours(1,:))
%         
%         plot(centres,cur_mean_b(:,j)-cur_std_b(:,j),'Color',plot_colours(2,:),'LineWidth',1)
%         plot(centres,cur_mean_b(:,j)+cur_std_b(:,j),'Color',plot_colours(2,:),'LineWidth',1)
%         % xx = 1:100;
%         h = fill([centres,fliplr(centres)],[(cur_mean_b(:,j)-cur_std_b(:,j))',fliplr((cur_mean_b(:,j)+cur_std_b(:,j))')],plot_colours(2,:));
%         set(h,'facealpha',.1)
%         plot(centres,cur_mean_b(:,j),'LineWidth',2,'Color',plot_colours(2,:))
%         axis square;
%     end
%     
% end

boot_samps = 100;
sub_sample = false;
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

% figure
for i = 1:2
    for j = 1:2        
        subplot(2,2,(i-1)*2+j)
        % plot(centres,squeeze(CIs_all(1,j,:,i)),'Color',plot_colours(1,:),'LineWidth',1)
        hold on
        % plot(centres,squeeze(CIs_all(2,j,:,i)),'Color',plot_colours(1,:),'LineWidth',1)
        h = fill([centres,fliplr(centres)],[squeeze(CIs_all(1,j,:,i))',fliplr(squeeze(CIs_all(2,j,:,i))')],plot_colours(1,:),'LineStyle','none');
        set(h,'facealpha',.3)
        
        plot(centres,squeeze(x_binned_means(:,i,j)),'Color',plot_colours(1,:),'LineWidth',2)
        
        % plot(centres,squeeze(CIs_all(1,j+2,:,i)),'Color',plot_colours(2,:),'LineWidth',1)
        hold on
        % plot(centres,squeeze(CIs_all(2,j+2,:,i)),'Color',plot_colours(2,:),'LineWidth',1)
        h = fill([centres,fliplr(centres)],[squeeze(CIs_all(1,j+2,:,i))',fliplr(squeeze(CIs_all(2,j+2,:,i))')],plot_colours(2,:),'LineStyle','none');
        set(h,'facealpha',.3)
        
        plot(centres,squeeze(x_binned_means(:,i,j+2)),'Color',plot_colours(2,:),'LineWidth',2)
        
        if i == 2
            % ylim([-2,2]);
            yline(0,'--','LineWidth',2);
        else
            % insert ylim for pitch
            % ylim([0,40]);
        end
        box off
    end
end

ax11 = subplot(2,2,1);
title("Left trials")
% xticklabels([])
% xlabel("Linearized position (cm)")
ylabel(["Forward velocity"; "(cm/s)"])
% ylim([0,40])
% axis('square')

ax12 = subplot(2,2,2);
title("Right trials")
% xticklabels([])
% xlabel("Linearized Position (cm)")
% ylim([0,40])
yticklabels([])
% axis('square')

ax13 = subplot(2,2,3);
%title("Left Trials")
ylabel(["Angular velocity"; "(rad/s)"])
xlabel("Linearized position (cm)")
% ylim([-2,2])
% axis('square')

ax14 = subplot(2,2,4);
%title("Right Trials")
xlabel("Linearized position (cm)")
% ylim([-2,2])
yticks([-0.5,0,0.5,1])
yticklabels([])
% axis('square')

linkaxes([ax11,ax12])
linkaxes([ax13,ax14])


%% Size and positioning

