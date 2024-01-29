function [dot_prods,dot_prods_yaw,h_boots] = plot_va_yaw_weights_compare(Wout_online_cell,bmi_weights_cell)
% 29/06/2023

% Plot example weights, and results of dot products of weights.
% Possibly be good to also add the weights only in neuron masks.
num_mice = size(bmi_weights_cell,1);
num_days = size(bmi_weights_cell,2);

ex_md = [5,1];
yaw_ind = 4;
nmodels = 5;

% Plot image of weights
ex_w = squeeze(Wout_online_cell{ex_md(1)}(:,2));
exbmi_w = squeeze(bmi_weights_cell{ex_md(1),ex_md(2)}(1,:,yaw_ind));

figure
% Transpose needed for mice 1-3
cur_weights = reshape(ex_w(1:end-1),128,128);

up_lim = max(cur_weights(:));
low_lim = min(cur_weights(:));

lim_val = max(abs([low_lim,up_lim]));

% [pix_neuron_list,pix_neuron_mat,sorted_pix_list,sorted_pix_is_neuron,ex_mask] = assign_pixel_weights_to_neurons(cur_weights,ex_stat);
% [pix_neuron_list,pix_neuron_mat,sorted_pix_list,sorted_pix_is_neuron,ex_mask] = assign_pixel_weights_to_neurons_11042023(cur_weights,ex_stat);
% Weights
% ax1 = subplot(plot_rows,4,6+4);
imagesc(cur_weights,[-lim_val,lim_val])
axis('square')
box on
% axis equal
% axis off
xticks([])
yticks([])
colormap('redblue')
colorbar
title("View Angle Weights")

figure
% Transpose needed for mice 1-3. Not transpose for 4-5
% Invert as beta is negative
cur_weights = -1.0*reshape(exbmi_w(1:end-1),128,128);
% cur_weights = reshape(exbmi_w(1:end-1),128,128);

up_lim = max(cur_weights(:));
low_lim = min(cur_weights(:));

lim_val = max(abs([low_lim,up_lim]));

% [pix_neuron_list,pix_neuron_mat,sorted_pix_list,sorted_pix_is_neuron,ex_mask] = assign_pixel_weights_to_neurons(cur_weights,ex_stat);
% [pix_neuron_list,pix_neuron_mat,sorted_pix_list,sorted_pix_is_neuron,ex_mask] = assign_pixel_weights_to_neurons_11042023(cur_weights,ex_stat);
% Weights
% ax1 = subplot(plot_rows,4,6+4);
imagesc(cur_weights,[-lim_val,lim_val])
axis('square')
box on
% axis equal
% axis off
xticks([])
yticks([])
colormap('redblue')
colorbar
title("Ball Angular Velocity Weights")

%% Calculate dot product and plot
% NB: Do I need to transpose the online weights for m1-3 but not for m4-5
% here?? ARRGGH

dot_prods = nan.*ones(nmodels,num_mice,num_days);
dot_prods_yaw = nan.*ones(nmodels,nmodels,num_mice,num_days);

for m = 1:num_mice
    Wout  = squeeze(Wout_online_cell{m}(1:end-1,2));
    % Add in transpose for m1-3
    if m < 4
        ww = reshape(Wout,128,128)';
        Wout = ww(:);
    end
    for d = 1:num_days
        if ~isempty(bmi_weights_cell{m,d})
            cur_b = bmi_weights_cell{m,d};
            for n = 1:nmodels
                % cur_n = squeeze(cur_b(n,1:end-1,yaw_ind));
                cur_n = -1.0*squeeze(cur_b(n,1:end-1,yaw_ind));
                dot_prods(n,m,d) = sum((cur_n.*Wout'))/(norm(cur_n)*norm(Wout));
                for i = 1:nmodels
                    if i ~= n 
                        cur_ni = -1.0*squeeze(cur_b(i,1:end-1,yaw_ind));
                        dot_prods_yaw(n,i,m,d) = sum((cur_n.*cur_ni))/(norm(cur_n)*norm(cur_ni));
                    end
                end
            end
                
            
        end
    end
end

%% Plot
num_sess = num_mice*num_days;
plot_off = linspace(-0.4,0.4,num_sess);

figure
mean_dot_prods = squeeze(mean(dot_prods,'omitnan'));
scatter(plot_off+ones(1,num_sess),mean_dot_prods(:)','filled','k')
hold on
plot([plot_off(1)+1,plot_off(end)+1],[mean(mean_dot_prods(:),'omitnan'),mean(mean_dot_prods(:),'omitnan')],'k','LineWidth',2)

title(["View Angle and Ball Angular Velocity"; "Weight Comparison"])
ylabel("Cosine Similarity")
yline(0,'--','LineWidth',2);
xticks([1])
xticklabels(["View Angle Vs Ball Angular Velocity"])
ylim([-0.5,0.5])
xlim([0.5,1.5])

% yaw with other yaw same session
figure
mean_yaw_prod = squeeze(mean(dot_prods_yaw,[1,2],'omitnan'));
scatter(plot_off+ones(1,num_sess),squeeze(mean_yaw_prod(:))','filled','k')
hold on
plot([plot_off(1)+1,plot_off(end)+1],[mean(squeeze(mean_yaw_prod(:)),'omitnan'),mean(squeeze(mean_yaw_prod(:)),'omitnan')],'k','LineWidth',2)

title(["View Angle and Ball Angular Velocity"; "Weight Comparison"])
ylabel("Cosine Similarity")
yline(0,'--','LineWidth',2);
xticks([])

figure
mean_dot_prods = squeeze(mean(dot_prods,'omitnan'));
scatter(plot_off+ones(1,num_sess),mean_dot_prods(:)','filled','k')
hold on
plot([plot_off(1)+1,plot_off(end)+1],[mean(mean_dot_prods(:),'omitnan'),mean(mean_dot_prods(:),'omitnan')],'k','LineWidth',2)

scatter(plot_off+2.*ones(1,num_sess),squeeze(mean_yaw_prod(:))','filled','k')
hold on
plot([plot_off(1)+2,plot_off(end)+2],[mean(squeeze(mean_yaw_prod(:)),'omitnan'),mean(squeeze(mean_yaw_prod(:)),'omitnan')],'k','LineWidth',2)

title(["View Angle and Ball Angular Velocity"; "Weight Comparison"])
ylabel("Cosine Similarity")
yline(0,'--','LineWidth',2);
xticks([1,2])
% xlabel("Weights Compared")
xticklabels([])
axis('square')

% Get mean and sem
dots_ready = nan.*ones(2,num_mice,num_days);
for m = 1:num_mice
    cur_dots = mean_dot_prods(m,:);
    cur_dots_yaw = mean_yaw_prod(m,:);
    dots_ready(1,m,1:sum(~isnan(cur_dots))) = cur_dots(~isnan(cur_dots));
    dots_ready(2,m,1:sum(~isnan(cur_dots_yaw))) = cur_dots_yaw(~isnan(cur_dots_yaw));
end

% boot_samps = 1000;
% num_trials = 4;
% 
% p_boots = NaN(2,1);
% 
% [p_boots(1), bootstats, bootstats_center, bootstats_sem] = get_bootstrap_results_equalsamples(squeeze(dots_ready(1,:,:)),squeeze(dots_ready(2,:,:)),boot_samps,num_trials,'mean');

[all_p_boot,all_centres,all_sems] = run_H_boot_ets(squeeze(dots_ready(1,:,:)), squeeze(dots_ready(2,:,:)),false);

h_boots.all_p_boot = all_p_boot;
h_boots.all_centres = all_centres;
h_boots.all_sems = all_sems;
