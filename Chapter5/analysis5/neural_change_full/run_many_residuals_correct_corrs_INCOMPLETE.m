function [] = run_many_residuals_correct_corrs_INCOMPLETE(z_binned_res_cell,z_binned_cell,x_binned_cell,error_cell,tbt_cell,n_W_online_cell,centres,error_thresh)
% 12/09/2023
% UNFINISHED: SWITCHED TO USING SAMPLES
num_mice = size(x_binned_cell,1);
num_days = size(x_binned_cell,2);

types_vec = [1,4,7,10];
nbins = 50;
% First project binned residuals through decoder weights

% Weights cells already weights in neurons from pixels.

projected_cell = cell(num_mice,num_days);

mean_binned = NaN(nbins,num_mice,num_days);
mean_all = NaN(num_mice,num_days);

mean_binned_b = NaN(nbins,num_mice,num_days);
mean_all_b = NaN(num_mice,num_days);

ro_all = nan.*ones(num_mice,num_days);
ro_bins = nan.*ones(nbins,num_mice,num_days);
ro_bins_shuff = nan.*ones(num_shuffles,nbins,num_mice,num_days);

for m = 1:num_mice
    
    for d = 1:num_days
        if ~isempty(x_binned_cell{m,d})
            Wout = n_W_online_cell{m,d};
            tbt_details = tbt_cell{m,d};
            zz = z_binned_res_cell{m,d};
            zb = z_binned_cell{m,d};
            % Don't include the bias term
            on_output = NaN(size(zz,1),size(zz,2));
            on_output_b = NaN(size(zz,1),size(zz,2));
            for t = 1:size(zz,1)
                on_output(t,:) = squeeze(zz(t,:,:))*Wout;
                on_output_b(t,:) = squeeze(zb(t,:,:))*Wout;
            end
            projected_cell{m,d} = on_output;
            
            cur_trials = find(ismember(tbt_details(3,:),types_vec([3,4])));
            
            % Calculate mean abs projected output
            mean_binned(:,m,d) = mean(abs(on_output(cur_trials,:)));
            mean_all(m,d) = mean(abs(on_output(cur_trials,:)),'all');
            
            mean_binned_b(:,m,d) = mean(abs(on_output_b(ismember(tbt_details(3,:),types_vec([3,4])),:)));
            mean_all_b(m,d) = mean(abs(on_output_b(ismember(tbt_details(3,:),types_vec([3,4])),:)),'all');
            
            
            % get correlations with corrective movements
            error_mat = error_cell{m,d};
            x_binned = x_binned_cell{m,d};
            % tbt_details = tbt_cell{m,d};

            
            
            cur_errors = error_mat(cur_trials,:);
            % cur_yaw = squeeze(x_binned(cur_trials,:,5)); % Yaw from decoded view angle
            % TODO: Calculate decoder yaw from binned decoder view angle.
            
            ro = corr(cur_errors(abs(cur_errors)>error_thresh),cur_yaw(abs(cur_errors)>error_thresh)); 
            ro_all(m,d) = ro;
            for b = 1:nbins
                ro_bins(b,m,d) = corr(cur_errors(abs(cur_errors(:,b))>error_thresh,b),cur_yaw(abs(cur_errors(:,b))>error_thresh,b));
            end
            
            % calculate shuffled corelations
            for n = 1:num_shuffles
                cur_yaw_shuff = cur_yaw(:);
                cur_yaw_shuff = cur_yaw_shuff(randperm(length(cur_yaw_shuff)));
                cur_yaw_shuff = reshape(cur_yaw_shuff,size(cur_yaw,1),size(cur_yaw,2));
                for b = 1:nbins
                    ro_bins_shuff(n,b,m,d) = corr(cur_errors(abs(cur_errors(:,b))>error_thresh,b),cur_yaw_shuff(abs(cur_errors(:,b))>error_thresh,b));
                end
            end
             
        end
    end
end
