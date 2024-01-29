function [sig_results_all,sig_check_all] = run_many_sig_tuning_check(z_cell, virmen_cell, tbt_cell, num_shuffles, linearise_x, nbins, boot_samps, shuff_limit, xnum)
% 29/03/2022

% 

num_mice = size(z_cell,1);
num_days = size(z_cell,2);
sig_results_all = cell(num_mice,num_days);
sig_check_all = cell(num_mice,num_days);

for m = 1:num_mice
    for d = 1:num_days
        if ~isempty(z_cell{m,d})
            % [bootstrap_means,sig_results,sig_check,centres] = significant_tuning_check_v2(z_cell{m,d}, virmen_cell{m,d}, tbt_cell{m,d}, num_shuffles, linearise_x, nbins, boot_samps, shuff_limit, xnum);
            [bootstrap_means,sig_results,sig_check,centres] = significant_tuning_check_v3_24102022(z_cell{m,d}, virmen_cell{m,d}, tbt_cell{m,d}, num_shuffles, linearise_x, nbins, boot_samps, shuff_limit, xnum);
            sig_results_all{m,d} = sig_results;
            sig_check_all{m,d} = sig_check;
        end
    end
end