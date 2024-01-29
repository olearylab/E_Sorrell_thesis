function [bootstats_center,bootstats_sem] = plot_many_KH_corrs(full_pairs_means_cell)
% 05/07/2023

% Plot example session clustered neurons
% Plot summary data...

ex_md = [4,3];
ex_pairs = full_pairs_means_cell{ex_md(1),ex_md(2)};
titles = ["Ball Left";"Ball Right";"BMI Left";"BMI Right"];
T = clusterdata(squeeze(ex_pairs(1,:,:)),1);
for i = 1:4
    % T = clusterdata(squeeze(ex_pairs(i,:,:)),1);
    [a,inds] = sort(T);
    figure 
    imagesc(squeeze(ex_pairs(i,inds,inds)),[-1,1])
    colormap('redblue')
    xticks([])
    yticks([])
    title(titles(i));
    axis('square')
end

% compare sorting
T = clusterdata(squeeze(ex_pairs(1,:,:)),1);
[a,inds] = sort(T);
figure 
imagesc(squeeze(ex_pairs(1,inds,inds)),[-1,1])
colormap('redblue')
xticks([])
yticks([])
figure 
imagesc(squeeze(ex_pairs(3,inds,inds)),[-1,1])
colormap('redblue')
xticks([])
yticks([])
% title(titles(i));
T = clusterdata(squeeze(ex_pairs(3,:,:)),1);
[a,inds] = sort(T);
figure 
imagesc(squeeze(ex_pairs(1,inds,inds)),[-1,1])
colormap('redblue')
xticks([])
yticks([])
figure 
imagesc(squeeze(ex_pairs(3,inds,inds)),[-1,1])
colormap('redblue')
xticks([])
yticks([])
%% Summary

num_mice = size(full_pairs_means_cell,1);
num_days = size(full_pairs_means_cell,2);

full_means = nan.*ones(4,num_mice,num_days);

full_sim = nan.*ones(4,4,num_mice,num_days);

for m = 1:num_mice
    for d = 1:num_days
        if ~isempty(full_pairs_means_cell{m,d})
            cur_pairs = full_pairs_means_cell{m,d};
            for i = 1:4
                cur_block = squeeze(cur_pairs(i,:,:));
        
                low_t = tril(nan.*ones(size(cur_block)));
                
                cur_block = cur_block + low_t;
                
                cur_block_ready = cur_block(:);
                full_means(i,m,d) = mean(cur_block_ready,'omitnan');
                for j = 1:4
                    cur_block_j = squeeze(cur_pairs(j,:,:));
        
                    low_t = tril(nan.*ones(size(cur_block_j)));
                
                    cur_block_j = cur_block_j + low_t;
                
                    cur_block_ready_j = cur_block_j(:);
                    
                    full_sim(i,j,m,d) = sum(cur_block_ready(~isnan(cur_block_ready)).*cur_block_ready_j(~isnan(cur_block_ready_j)))/(norm(cur_block_ready(~isnan(cur_block_ready)))*norm(cur_block_ready_j(~isnan(cur_block_ready_j))));
                    
                end
            end
        end
    end
end

num_sess = num_mice*num_days;
figure
plot_off = linspace(-0.4,0.4,num_sess);

for i = 1:4
    scatter(plot_off+i.*ones(1,num_sess),squeeze(full_means(i,:))','filled','k')
    hold on
    plot([i+plot_off(1),i+plot_off(end)],[mean(squeeze(full_means(i,:)),'omitnan'),mean(squeeze(full_means(i,:)),'omitnan')],'k','LineWidth',2)
end
% title("Noise Correlations")
ylabel("Noise Correlations")
xticks([1,2,3,4])
xtickangle(30)
xticklabels(["Ball Left";"Ball Right";"BMI Left";"BMI Right"])


figure
plot_off = linspace(-0.4,0.4,num_sess);

i = 1;

lr_combined_sims = mean([squeeze(full_sim(1,3,:))';squeeze(full_sim(2,4,:))'],'omitnan');


% for i = 1:1
scatter(plot_off+i.*ones(1,num_sess),lr_combined_sims,'filled','k')
hold on
plot([i+plot_off(1),i+plot_off(end)],[mean(lr_combined_sims,'omitnan'),mean(lr_combined_sims,'omitnan')],'k','LineWidth',2)
% end
title(["Similarity of";"Noise Correlations"])
ylabel("Cosine Similarity")
% xticks([1,2,3,4])
% xtickangle(30)
% xticklabels(["Ball Left";"Ball Right";"BMI Left";"BMI Right"])
ylim([0,1])
xlim([0.5,1.5])
xticklabels([])



%% Hierarchical bootstrap

lr_combined_sims = (squeeze(full_sim(1,3,:,:))+squeeze(full_sim(2,4,:,:)))/2;

stats_data = nan.*ones(num_mice,num_days);
for m = 1:num_mice

    cur_data = lr_combined_sims(m,:);

    stats_data(m,1:sum(~isnan(cur_data))) = cur_data(~isnan(cur_data));
end

rng(1);
boot_samps = 1000;
num_trials = 4;
bootstats = get_bootstrapped_equalsamples(stats_data,boot_samps,num_trials,'mean');
%Get mean and SEM of bootstrapped samples:
bootstats_sem = std(bootstats);
bootstats_center = mean(bootstats);