function [all_corrs,corr_val,p_val,h_boots] = compare_pix_neu_many(mice,data_days,rep,training_trials,stat_cell_CNN_offline)
% 30/05/2023

% calculate correlation between neuron weights, and neurons in pixel
% weights for many offline sessions
res_path = "../shared_files/offline_results/";
res_path_n = "../shared_files/offline_neuron_results_CNN_16112022/";
% just do for rep 1 of each

xvars = [1,7,3,4];

num_mice = length(mice);

num_mice = length(mice);
n_days = zeros(num_mice,1);
for m = 1:num_mice
    n_days(m) = length(data_days{m});
end
num_days = max(n_days);

all_corrs_mat = nan.*ones(num_mice,num_days,4);

all_corrs = nan.*ones(length(stat_cell_CNN_offline),length(xvars));
ind = 1;
for m = 1:num_mice
    data_d = data_days{m};
    num_days = length(data_d);
    for d = 1:num_days
        stat = stat_cell_CNN_offline{ind};
        results_struct = importdata(res_path + "res_struct_" + mice(m) + "_" + data_d(d) + "_rep"+ rep + "_" + training_trials + ".mat");
        Wout = results_struct.Wout(:,xvars);
        
        if m<4
            for i = 1:length(xvars)
                ww = reshape(Wout(1:end-1,i),128,128);
                ww = ww';
                Wout(1:end-1,i) = ww(:);
            end
        end
        
        [n_weights_max_abs,n_weights_mean] = calc_weights_from_mask_11042023(Wout,stat);
        
        results_struct_n = importdata(res_path_n + "res_struct_" + mice(m) + "_" + data_d(d) + "_rep"+ rep + "_" + training_trials + ".mat");
        Wout_n = results_struct_n.Wout(1:end-1,xvars);
        
        % Plot example
        if ind == 1
            [corr_val,p_val] = compare_pixel_to_neuron_weights_30052023(Wout(:,1),stat,results_struct_n.Wout(:,1));
        end
        
        cur_corrs = nan.*ones(length(xvars),1);
        for i = 1:length(xvars)
            cur_corrs(i) = corr(Wout_n(:,i),n_weights_mean(:,i));
        end
        all_corrs(ind,:) = cur_corrs;
        ind = ind+1;
        
        all_corrs_mat(m,d,:) = cur_corrs;
    end
end


% Plot summary
figure
violin(all_corrs,'facecolor',[0.5,0.5,0.5],'medc',[]);
ylim([0.5,1])
xticks([1,2,3,4])
yticks([0.5,0.75,1])
title(["Correlation Between Neuron Weights";"In Pixel and Neuron Decoders"])
ylabel("Correlation")
xticklabels(["Y Position";"View Angle";"Forward Velocity";"Angular Velocity"])
xtickangle(30)
axis('square')

%% Hierarchical Bootstrap

% now for neuron pixel comparison
% all_p_boot = nan.*ones(4,1);
all_centres = nan.*ones(4,1);
all_sems = nan.*ones(4,1);

boot_samps = 1000;
num_trials = 4;

rng(1) % set random seed
for i = 1:4
    [bootstats] = get_bootstrapped_equalsamples(squeeze(all_corrs_mat(:,:,i)),boot_samps,num_trials,'mean');
    all_sems(i) = std(bootstats);
    all_centres(i) = mean(bootstats);
end

% h_boots.all_p_boot = all_p_boot;
h_boots.all_centres = all_centres;
h_boots.all_sems = all_sems;
