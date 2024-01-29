function [z_binned_res_cell,z_binned_cell,x_binned_cell,full_all_means] = get_z_binned_res_cell_samps(z_cell,virmen_cell,tbt_cell,n_Wout_online_cell,model_params,normalise_z,nbins)
% 12/09/2023
% For online view angle weights
% 
num_mice = size(z_cell,1);
num_days = size(z_cell,2);

z_binned_res_cell = cell(num_mice,num_days);
z_binned_cell = cell(num_mice,num_days);
x_binned_cell = cell(num_mice,num_days);

full_all_means = NaN(2,num_mice,num_days);

for m = 1:num_mice
    for d = 1:num_days
        if ~isempty(z_cell{m,d})
            % x_vec = [6,7,15,17,on_ind,res_ind,vav_ind,vav_res_ind];
            % [z_binned] = calc_binned_residual_activity(z_cell{m,d},virmen_cell{m,d},tbt_cell{m,d},model_params,normalise_z,nbins);
            [z_binned,z_binned_res,x_binned,all_means] = calc_binned_residual_activity_from_samps(z_cell{m,d},virmen_cell{m,d},tbt_cell{m,d},n_Wout_online_cell{m,d},model_params,normalise_z,nbins,[],[]);

            z_binned_res_cell{m,d} = z_binned_res;
            z_binned_cell{m,d} = z_binned;
            x_binned_cell{m,d} = x_binned;
            full_all_means(:,m,d) = all_means;
                
        end
    end
end
