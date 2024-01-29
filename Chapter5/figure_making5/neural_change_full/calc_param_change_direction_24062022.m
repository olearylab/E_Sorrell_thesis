function [final_averages] = calc_param_change_direction_24062022(all_mean_curve_params,all_full_results_both,sig_results_all,sig_summary_cell)
% 24/06/2022

% Calculate whether changes are increases or decreases.
% Also average across direction, mice and session.

% Issue of potentially not adding up to same overall average.

% all_mean_curve_params can be all_full_param_means - this is better

% calculate proportion of neurons tuned to both, that have changes

%% 

num_mice = size(all_full_results_both,1);
num_days = size(all_full_results_both,2);

% Only keep significant neurons, with significant changes in parameter
all_sig_changes = cell(num_mice,num_days,2,3);
all_sig_num_changes = zeros(num_mice,num_days,2,2,3);
all_sig_proportions = nan.*ones(num_mice,num_days,2,2,3);

for m = 1:num_mice
    for d = 1:num_days
        if ~isempty(all_full_results_both{m,d})
            full_results = all_full_results_both{m,d};
            cur_sig = sig_results_all{m,d};
            cur_sig_summary = sig_summary_cell{m,d};
            mean_curve_params = all_mean_curve_params{m,d};
            for i = 1:2
                % cur_neurons = cur_sig(i,:)==1 && cur_sig(i+2,:)==1;
                cur_neurons = squeeze(cur_sig_summary(i,1,:));
                for j = 1:3
                    cur_param_neurons = squeeze(full_results(:,i,j));
                    all_sig_changes{m,d,i,j} = squeeze(mean_curve_params(i+2,cur_neurons&cur_param_neurons,j) - mean_curve_params(i,cur_neurons&cur_param_neurons,j));
                    all_sig_num_changes(m,d,i,1,j) = sum(all_sig_changes{m,d,i,j}>0);
                    all_sig_num_changes(m,d,i,2,j) = sum(all_sig_changes{m,d,i,j}<0);
                    tot_nums = squeeze(all_sig_num_changes(m,d,i,1,j)) + squeeze(all_sig_num_changes(m,d,i,2,j));
                    all_sig_proportions(m,d,i,:,j) = all_sig_num_changes(m,d,i,:,j)./tot_nums;
                    
                end
            end
        end
    end
end

% Could take average from sums (so each neuron weighted equally), or weight
% each session equally. Do the latter for now.
% give single number for each parameter which is proportion of changed
% neurons that have increased.
final_averages = squeeze(mean(all_sig_proportions(:,:,:,1,:),[1,2,3],'omitnan'));