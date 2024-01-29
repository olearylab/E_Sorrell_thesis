function [z_binned_res_cell,z_binned_cell,x_binned_cell,full_all_means] = get_z_binned_res_cell_samps_subs(z_cell,virmen_cell,tbt_cell,n_Wout_cell,model_params,normalise_z,nbins)
% 13/09/2023
% For many subsamples in training angular velocity weights
% Need to edit once I have to full results
% 
num_mice = size(z_cell,1);
num_days = size(z_cell,2);

% nsubs = size(n_Wout_cell{1,1},1);
% REMOVE
nsubs = 1;
nmodels = size(n_Wout_cell{1,1},2);
tot_reps = nsubs*nmodels;

z_binned_res_cell = cell(num_mice,num_days);
z_binned_cell = cell(num_mice,num_days);
x_binned_cell = cell(tot_reps,num_mice,num_days);

full_all_means = NaN(tot_reps,2,num_mice,num_days);

for m = 1:num_mice
    for d = 1:num_days
        if ~isempty(z_cell{m,d})
            % x_vec = [6,7,15,17,on_ind,res_ind,vav_ind,vav_res_ind];
            % [z_binned] = calc_binned_residual_activity(z_cell{m,d},virmen_cell{m,d},tbt_cell{m,d},model_params,normalise_z,nbins);
            
            cur_w = n_Wout_cell{m,d};
            cur_w = permute(cur_w,[3,1,2]);
            
            % REMOVE
            cur_w = squeeze(cur_w(:,30,:));
            % uncomment
            % cur_w = cur_w(:,:);
            
            for n = 1:tot_reps
                [z_binned,z_binned_res,x_binned,all_means] = calc_binned_residual_activity_from_samps(z_cell{m,d},virmen_cell{m,d},tbt_cell{m,d},cur_w(:,n),model_params,normalise_z,nbins,[],[]);
                if n == 1
                    z_binned_res_cell{m,d} = z_binned_res;
                    z_binned_cell{m,d} = z_binned;
                end
                x_binned_cell{n,m,d} = x_binned;
                full_all_means(n,:,m,d) = all_means;
            end
        end
    end
end
