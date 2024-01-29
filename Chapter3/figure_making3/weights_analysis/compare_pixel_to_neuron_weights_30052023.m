function [corr_val,p_val] = compare_pixel_to_neuron_weights_30052023(Wout,stat,Wout_CNN)
% 22/03/2023

% function for checking whether neurons are weighted in similar rank order
% when just neurons used or all pixels used.

% Potential issue of initialisation for first 3 mice.

% offline_index = [3,7];


[n_weights_max_abs,n_weights_mean] = calc_weights_from_mask_11042023(Wout,stat);
% n_weights_1 = n_weights_max_abs;
% n_weights_2 = abs(n_weights_mean);

% Sort neurons according to mean weight size
[sorted_weights_mean,I_weights_mean] = sort(abs(n_weights_mean));
% Sort neurons according to max weight size
[sorted_weights_max_abs,I_weights_max_abs] = sort(n_weights_max_abs);

[~,rank_weights_mean] = sort(I_weights_mean);
[~,rank_weights_max_abs] = sort(I_weights_max_abs);

% if ~isempty(zdata)
%     [n_weights_sum,n_weights_scaled_sum] = calc_pix_neuron_weight_importance_22032023(Wout,stat,zdata,mouse);
%     n_weights_1 = abs(n_weights_sum);
%     n_weights_2 = abs(n_weights_scaled_sum);
% end
% Compare max rank order to mean rank order

% % Sort neurons according to mean weight size
% [sorted_weights_1,I_weights_1] = sort(n_weights_1);
% % Sort neurons according to max weight size
% [sorted_weights_2,I_weights_2] = sort(n_weights_2);
% 
% [~,rank_weights_1] = sort(I_weights_1);
% [~,rank_weights_2] = sort(I_weights_2);

% figure
% for i = 1:size(rank_weights_1,2)
%     subplot(1,2,i)
%     scatter(rank_weights_1(:,i),rank_weights_2(:,i));
% end
% 
% figure
% for i = 1:size(rank_weights_1,2)
%     subplot(1,2,i)
%     scatter(rank_weights_2(:,i),rank_weights_mean(:,i));
% end
% 
% figure
% for i = 1:size(rank_weights_1,2)
%     subplot(1,2,i)
%     scatter(rank_weights_2(:,i),rank_weights_max_abs(:,i));
% end
% 
% weight_order = flipud(I_weights_2);

Wout_n = Wout_CNN(1:end-1);
% 
% if ~isempty(zdata_n)
%     initialise_params;
%     model_params.reg_images = false;
%     model_params.spatial = false;
%     create_reg = false;
%     disp("Pre-Processing")
%     [zdatafilt, train_mean, Rfixed] = preprocess_images_05_2021(zdata_n,model_params,[],[],[],create_reg);
%     disp("Pre-Processing Complete") 
%     z_stds = std(zdatafilt);
%     Wout_n = Wout_n.*z_stds';
% end
% 
[sorted_weights,I_weights] = sort(abs(Wout_n));
% 
[~,rank_weights] = sort(I_weights);
% 
% figure
% scatter(rank_weights,rank_weights_mean);

figure
scatter(Wout_n,n_weights_mean,'filled','MarkerFaceColor',[0.5,0.5,0.5],'MarkerEdgeColor',[0.5,0.5,0.5]);
hold on
xline(0,'--','LineWidth',2);
yline(0,'--','LineWidth',2);
axis('square')

[corr_val,p_val] = corr(Wout_n,n_weights_mean);

mdl = fitlm(Wout_n,n_weights_mean);
xx = linspace(min(Wout_n),max(Wout_n),100)';
[ypred,yci] = predict(mdl,xx);
plot(xx,yci(:,1),'--','LineWidth',2,'Color','k')
plot(xx,yci(:,2),'--','LineWidth',2,'Color','k')
plot(xx,ypred,'-','LineWidth',2,'Color','k')
h = fill([xx',fliplr(xx')],[yci(:,2)',fliplr(yci(:,1)')],'k','LineStyle','none');
set(h,'facealpha',.1)