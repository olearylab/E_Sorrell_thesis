function [] = plot_example_weights(ex_weights,ex_stat,ex_im,trans_weights)

%% Temporary alternate plot
plot_nums = [1,7,3,4]; % ypos, VA, pitch, yaw
ex_w_num = 1;
figure;
set(gcf,'position',[1,42,1440,254])
% Mean Downsampled Image (Maybe use full image?)
% ax = subplot(plot_rows,4,5+4);
ax = subplot(1,4,1);
mm = mean(ex_im(:));
ss = std(ex_im(:));
ll = mm - ss;
uu = mm + 5*ss;
imagesc(ex_im,[ll,uu])
axis('square')
box on
% axis equal
% axis off
xticks([])
yticks([])
colormap(ax,'gray')
title("Mean downsampled image")
set(gca,'Position',[0.13,0.11,0.126394067899881,0.815])
hold on
% scalebar (1.2 microns per pixel)
plot([10,10+(100/4.8)],[15,15],'LineWidth',4,'Color','w')
% colorbar
   
if trans_weights
    cur_weights = reshape(ex_weights(1:end-1,plot_nums(ex_w_num)),128,128)';
else
    cur_weights = reshape(ex_weights(1:end-1,plot_nums(ex_w_num)),128,128);
end

up_lim = max(cur_weights(:));
low_lim = min(cur_weights(:));

lim_val = max(abs([low_lim,up_lim]));

[pix_neuron_list,pix_neuron_mat,sorted_pix_list,sorted_pix_is_neuron,ex_mask] = assign_pixel_weights_to_neurons_11042023(cur_weights,ex_stat);
% Weights
% ax1 = subplot(plot_rows,4,6+4);
ax1 = subplot(1,4,2);
imagesc(cur_weights,[-lim_val,lim_val])
axis('square')
box on
% axis equal
% axis off
xticks([])
yticks([])
colormap(ax1,'redblue')
colorbar
title("Y position weights")

cur_n_weights = cur_weights;
cur_n_weights(ex_mask==0) = nan;
% Only weights within neurons
% ax2 = subplot(plot_rows,4,7+4);
ax2 = subplot(1,4,3);
h = imagesc(cur_n_weights,[-lim_val,lim_val]);
set(h, 'AlphaData', ~isnan(cur_n_weights))
% axis off
axis('square')
box on
% axis equal
xticks([])
yticks([])
colormap(ax2,'redblue')
set(gca,'Color',[0.8,0.8,0.8])
colorbar
title("Weights in neurons")


cur_o_weights = cur_weights;
cur_o_weights(ex_mask==1) = nan;

% ax3 = subplot(plot_rows,4,8+4);
ax3 = subplot(1,4,4);
h = imagesc(cur_o_weights,[-lim_val,lim_val]);
set(h, 'AlphaData', ~isnan(cur_o_weights))
% axis off
axis('square')
box on
% axis equal
xticks([])
yticks([])
colormap(ax3,'redblue')
set(gca,'Color',[0.8,0.8,0.8])
colorbar
title("Weights outside neurons")


%% Plot zoomed in version
% Just adjust axes limits of main figure
figure
set(gcf,'position',[1,42,1440,254])
% Mean Downsampled Image (Maybe use full image?)
% ax = subplot(plot_rows,4,5+4);
ax = subplot(1,4,1);
mm = mean(ex_im(:));
ss = std(ex_im(:));
ll = mm - ss;
uu = mm + 5*ss;
imagesc(ex_im,[ll,uu])
axis('square')
box on
% axis equal
% axis off
xticks([])
yticks([])
colormap(ax,'gray')
title("Mean downsampled image")
set(gca,'Position',[0.13,0.11,0.126394067899881,0.815])
hold on
% scalebar (1.2 microns per pixel)
plot([10,10+(50/4.8)],[15,15],'LineWidth',4,'Color','w')

% colorbar
ax1 = subplot(1,4,2);
imagesc(cur_weights,[-lim_val,lim_val])
axis('square')
box on
% axis equal
% axis off
xticks([])
yticks([])
colormap(ax1,'redblue')
colorbar
title("Y position weights")

cur_n_weights = cur_weights;
cur_n_weights(ex_mask==0) = nan;
% Only weights within neurons
% ax2 = subplot(plot_rows,4,7+4);
ax2 = subplot(1,4,3);
h = imagesc(cur_n_weights,[-lim_val,lim_val]);
set(h, 'AlphaData', ~isnan(cur_n_weights))
% axis off
axis('square')
box on
% axis equal
xticks([])
yticks([])
colormap(ax2,'redblue')
set(gca,'Color',[0.8,0.8,0.8])
colorbar
title("Weights in neurons")


cur_o_weights = cur_weights;
cur_o_weights(ex_mask==1) = nan;

% ax3 = subplot(plot_rows,4,8+4);
ax3 = subplot(1,4,4);
h = imagesc(cur_o_weights,[-lim_val,lim_val]);
set(h, 'AlphaData', ~isnan(cur_o_weights))
% axis off
axis('square')
box on
% axis equal
xticks([])
yticks([])
colormap(ax3,'redblue')
set(gca,'Color',[0.8,0.8,0.8])
colorbar
title("Weights outside neurons")

linkaxes([ax,ax1,ax2,ax3])
xlim([0.5,64.5])
ylim([0.5,64.5])