function [] = plot_example_curves_06122022(bootstrap_means,centres,example_neurons,ex_neu_dir,CI_vals)
% 20/09/2022

% Function for plotting some example tuning curves (for poster/paper)

% colours_vec = [[0 0 0]; [0 0.4470 0.7410]];
colours_vec = [[0.4940 0.1840 0.5560];[0.4660 0.6740 0.1880]];
num_examples = length(example_neurons);
zdim = size(bootstrap_means,4);
t_types = [1,4,7,10];
%%

z_binned_means = zeros(size(bootstrap_means,3),zdim,length(t_types));
  
for i = 1:length(t_types)
    z_binned_means(:,:,i) = mean(squeeze(bootstrap_means(i,:,:,:)),'omitnan');
end

% CIs, upper/lower x ntypes x nbins x zdim
CIs_all = zeros(2,length(t_types),size(bootstrap_means,3),zdim);
for i = 1:size(bootstrap_means,3)
    for j = 1:length(t_types)
        CIs_all(:,j,i,:) = prctile(squeeze(bootstrap_means(j,:,i,:)),CI_vals);
    end
end

%% Plot 
figure

for i = 1:num_examples

    subplot(ceil(num_examples/2),2,i)

    t_ind = ex_neu_dir(i);
    cur_n = example_neurons(i);
        
    for t = 1:2
        % error bars
        % plot(centres,z_binned_means(:,(p-1)*num_neurons+i,t)+squeeze(std(bootstrap_means(t,:,:,(p-1)*num_neurons+i))),'Linewidth',1,'Color',colours_vec(t,:))
        % plot(centres,squeeze(CIs_all(1,t_ind + 2*(t-1),:,cur_n)),'Linewidth',1,'Color',colours_vec(t,:))
        hold on
        % plot(centres,squeeze(CIs_all(2,t_ind + 2*(t-1),:,cur_n)),'Linewidth',1,'Color',colours_vec(t,:))
        % plot(centres,z_binned_means(:,(p-1)*num_neurons+i,t)-squeeze(std(bootstrap_means(t,:,:,(p-1)*num_neurons+i))),'Linewidth',1,'Color',colours_vec(t,:))
        h = fill([centres,fliplr(centres)],[squeeze(CIs_all(1,t_ind + 2*(t-1),:,cur_n))',fliplr(squeeze(CIs_all(2,t_ind + 2*(t-1),:,cur_n))')],colours_vec(t,:),'LineStyle','none');
        set(h,'facealpha',.3)
        % Means

        plot(centres,z_binned_means(:,cur_n,t_ind + 2*(t-1)),'Linewidth',3,'Color',colours_vec(t,:))

    end
    box off
    %yline(cur_sig_check(cur_n),'--','LineWidth',2);
    xticks([])
    %yticks([])
    % ylabel("Activity (a.u.)")
    axis('square')
end
subplot(ceil(num_examples/2),2,1)
ylabel(["Normalised"; "Activity (a.u.)"])
subplot(ceil(num_examples/2),2,3)
ylabel(["Normalised"; "Activity (a.u.)"])
subplot(ceil(num_examples/2),2,3)
xlabel("Linearised Position (cm)")
xticks([0,200])
subplot(ceil(num_examples/2),2,4)
xlabel("Linearised Position (cm)")
xticks([0,200])
