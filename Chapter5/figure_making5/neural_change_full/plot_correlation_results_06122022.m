function [] = plot_correlation_results_06122022(all_corrs_results,centres)
% 12/10/2022

% 

num_mice = size(all_corrs_results.all_nn_corrs,1);
num_days = size(all_corrs_results.all_nn_corrs,2);
nbins = size(all_corrs_results.all_binned_corrs{1,1},3);

% centres = centres*0.74;

colours_vec = [[0 0.4470 0.7410];[0.8500 0.3250 0.0980];[0.9290 0.6940 0.1250];[0.4940 0.1840 0.5560];[0.4660 0.6740 0.1880]];

tt_corrs_mat = nan.*ones(num_mice,num_days,4,4);
binned_corrs_mat = nan.*ones(num_mice,num_days,4,4,nbins);

for m = 1:num_mice
    for d = 1:num_days
        if ~isempty(all_corrs_results.all_nn_corrs{m,d})
            
            nn_corrs = all_corrs_results.all_nn_corrs{m,d};
            nn_noise_corrs = all_corrs_results.all_nn_noise_corrs{m,d};
            nn_res_corrs = all_corrs_results.all_nn_res_corrs{m,d};
            neuron_tt_corrs = all_corrs_results.all_neuron_tt_corrs{m,d};
            KH_binned_corrs = all_corrs_results.all_KH_binned_corrs{m,d};
            binned_corrs = all_corrs_results.all_binned_corrs{m,d};

            % Make figure like original correlation comparisons
            tt_corrs_mat(m,d,:,:) = mean(neuron_tt_corrs,3,'omitnan');
            binned_corrs_mat(m,d,:,:,:) = binned_corrs;
            
        end
    end
end

% %% Plot results
% % Single point per session, split mice horizontally slightly
% figure
% for m = 1:num_mice
%     for i = 1:4
%         scatter((i+0.1*m).*ones(num_days,1),tt_corrs_mat(m,:,i,i),'MarkerFaceColor',colours_vec(m,:),'MarkerEdgeColor',colours_vec(m,:))
%         hold on
%     end
%     ylim([0,1])
% end
% 
% 
% % Single point per session, split mice horizontally slightly.
% % Compare within to between sessions.
% figure
% for m = 1:num_mice
%     for i = 1:4
%         subplot(2,2,i)
%         for j = 1:4
%             scatter((j+0.1*m).*ones(num_days,1),tt_corrs_mat(m,:,i,j),'MarkerFaceColor',colours_vec(m,:),'MarkerEdgeColor',colours_vec(m,:))
%             hold on
%         end
%         ylim([0,1])
%     end
% end
% 
% % Example distribution across neurons
% % figure
% % for m = 1:num_mice
% %     for d = 1:num_days
% %         if ~isempty(all_corrs_results.all_nn_corrs{m,d})
% %             
% %             neuron_tt_corrs = all_corrs_results.all_neuron_tt_corrs{m,d};
% %             
% %             subplot(num_mice,num_days,(m-1)*num_days + d)
% %             boxplot(squeeze(neuron_tt_corrs(1,:,:))')
% %             
% %         end
% %     end
% % end
%   
% % Per bin correlations within type
plot_line = ["-";"--"];
figure
for m = 1:num_mice
    for i = 1:2
        subplot(1,2,i)
        for j = 1:2
            plot(centres,squeeze(mean(binned_corrs_mat(m,:,i+2*(j-1),i+2*(j-1),:),2,'omitnan')),'LineStyle',plot_line(j),'Color',colours_vec(m,:),'LineWidth',2)
            hold on
        end
        ylim([0,1])
    end
end
% 
% % Per bin correlations across type
% % Is there a specific location where differences are stronger or not?
% % Not bother with other across type checks
plot_line = ["-";"--"];
matched_index = [1,3;2,4;3,1;4,2];
figure
for m = 1:num_mice
    for i = 1:4
        subplot(2,2,i)
        for j = 1:2
            plot(centres,squeeze(mean(binned_corrs_mat(m,:,i,matched_index(i,j),:),2,'omitnan')),'LineStyle',plot_line(j),'Color',colours_vec(m,:),'LineWidth',2)
            hold on
        end
        ylim([0,1])
    end
end
% 
% % Plot per bin averaged across bins within types
% figure
% for m = 1:num_mice
%     for i = 1:4
%         scatter((i+0.1*m).*ones(num_days,1),squeeze(mean(binned_corrs_mat(m,:,i,i,:),5,'omitnan')),'MarkerFaceColor',colours_vec(m,:),'MarkerEdgeColor',colours_vec(m,:))
%         hold on
%     end
%     ylim([0,1])
% end
% 
% % Plot all type combinations
% x_vals = 1:4; % 
% % x_vals = x_vals + 0.1*median(1:num_mice);
% titles = ["Ball Left";"Ball Right";"BMI Left";"BMI Right"];
% figure
% for i = 1:4
%     subplot(2,2,i)
%     bar(x_vals,squeeze(mean(binned_corrs_mat(:,:,i,:,:),[1,2,5],'omitnan')),'w','LineWidth',2)
%     hold on
%     for m = 1:num_mice
%         for j = 1:4
%             scatter(((j+0.1*m)-0.3).*ones(num_days,1),squeeze(mean(binned_corrs_mat(m,:,i,j,:),5,'omitnan')),'MarkerFaceColor',colours_vec(m,:),'MarkerEdgeColor',colours_vec(m,:))
%         end
%     end
%     ylim([0,0.5])
%     title(titles(i))
%     xticks([1,2,3,4])
% end
% subplot(2,2,1)
% ylabel("Correlation")
% 
% subplot(2,2,3)
% ylabel("Correlation")

%% Correlation matrix 
% Only means in this case - population correlations
corr_image = squeeze(mean(binned_corrs_mat,[1,2,5],'omitnan'));
figure
imagesc(corr_image)
xticks([1,2,3,4])
yticks([1,2,3,4])
xticklabels(["Ball Left";"Ball Right";"BMI Left";"BMI Right"]);
yticklabels(["Ball Left";"Ball Right";"BMI Left";"BMI Right"]);
xtickangle(30);
title("Population Correlations")
colormap gray
axis('square')

%% trial-trial correlation matrix
corr_image = squeeze(mean(tt_corrs_mat,[1,2],'omitnan'));
figure
imagesc(corr_image)
xticks([1,2,3,4])
yticks([1,2,3,4])
xticklabels(["Ball Left";"Ball Right";"BMI Left";"BMI Right"]);
yticklabels(["Ball Left";"Ball Right";"BMI Left";"BMI Right"]);
xtickangle(30);
title("trial-trial Correlations")
colormap gray
axis('square')

%% Plot population correlations by bin
