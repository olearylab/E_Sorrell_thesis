function [all_means_cell,all_means_cell_ones,all_stds_cell] = run_many_decoder_dimension_check_n(n_W_online_cell,n_W_offline_cell,z_cell,virmen_cell,tbt_cell,model_params)
% 27/07/2023

% for weights in neurons

num_mice = size(virmen_cell,1);
num_days = size(virmen_cell,2);

all_means_cell = cell(num_mice,num_days);
all_stds_cell = cell(num_mice,num_days);

all_means_cell_ones = cell(num_mice,num_days);

for m = 1:num_mice
    for d = 1:num_days
        if ~isempty(virmen_cell{m,d})
            W_on = n_W_online_cell{m,d};
            W_off = n_W_offline_cell{m,d};
            disp("Running mouse " + m + " Day " + d)
            [all_means_cell{m,d},all_stds_cell{m,d}] = decoder_dimension_check_n(W_on,W_off,z_cell{m,d},virmen_cell{m,d},tbt_cell{m,d},model_params);
            
            % all pass weights check
            W_ones = ones(length(W_on),1);
            [all_means_cell_ones{m,d}] = decoder_dimension_check_n(W_ones,W_ones,z_cell{m,d},virmen_cell{m,d},tbt_cell{m,d},model_params);
            
        end
    end
end
