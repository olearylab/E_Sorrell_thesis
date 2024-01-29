function [thresh_num,thresh_num_RMSE,prop_90] = plot_decoder_dimensionality_and_n_weights_11042023(full_all_res,Wout,thresh_pct,stat)
% 14/10/2022

% Plot deocder performance by weights included, alongside weight size

[sorted_w,w_inds] = sort(abs(Wout(1:end-1)),'descend');

% Calculate some thresholds?
thresh_num = find(full_all_res(2,:)>thresh_pct*max(full_all_res(2,:)),1);
RMSE_range = max(full_all_res(1,:))-min(full_all_res(1,:));
thresh_num_RMSE = find(full_all_res(1,:)<max(full_all_res(1,:)) - thresh_pct*RMSE_range,1);

% RMSE
figure
% subplot(1,2,1)
yyaxis left
plot(full_all_res(1,:).*0.01,'LineWidth',2)

xline(thresh_num_RMSE,'--','LineWidth',2);
yline(full_all_res(1,thresh_num_RMSE).*0.01,'--','LineWidth',2);
ylabel("RMSE (m)")

xlim([0,size(full_all_res,2)])

yyaxis right
plot(sorted_w,'LineWidth',2)
ylabel("Weight")
xlabel("Number of Weights Included")
axis('square')

%%
% Weights must be transposed in order to match up with stat
W1 = reshape(Wout(1:end-1),128,128)';

[pix_neuron_list,pix_neuron_mat,sorted_pix_list,sorted_pix_is_neuron,final_footprints] = assign_pixel_weights_to_neurons_11042023(W1,stat);

[sorted_w,w_inds] = sort(abs(W1(:)),'descend');

%% Plot running proportion of weights that are in neruons
% Sorted from largest to smallest weight

running_prop = zeros(length(sorted_pix_is_neuron),1);

for i = 1:length(sorted_pix_is_neuron)
    
    running_prop(i) = sum(sorted_pix_is_neuron(1:i)==1)./i;
    
end

prop_90 = running_prop(thresh_num_RMSE);
figure
% subplot(1,2,2)
plot(running_prop,'k','LineWidth',2)
hold on
xline(thresh_num_RMSE,'--','LineWidth',2);
yline(prop_90,'--','LineWidth',2);
ylabel("Fraction of Weights in Neurons")
xlabel("Number of Weights Included")
ylim([0,1])
xlim([0,size(full_all_res,2)])
axis('square')

% % R^2
% figure
% yyaxis left
% plot(full_all_res(2,:),'LineWidth',2)
% 
% xline(thresh_num,'LineWidth',2);
% yline(full_all_res(2,thresh_num),'LineWidth',2);
% 
% yyaxis right
% plot([nan;sorted_w],'LineWidth',2)