function [va_corrs,yaw_corrs] = compare_weights_to_changes(Wout_online_cell,bmi_weights_cell,stat_cell,all_overall_n_means,z_stds)
% 26/07/2023

% Plot correlation between decoder weights and increases in neural activity

% Also compare changes to movement subspace (trained offline)
yaw_ind = 4;
W_ind = 2;

num_mice = size(bmi_weights_cell,1);
num_days = size(bmi_weights_cell,2);
nmodels = 5;

% 3 for raw, scaled, and scaling factor
va_corrs = nan.*ones(3,num_mice,num_days);
yaw_corrs = nan.*ones(3,nmodels,num_mice,num_days);

for m = 1:num_mice
    Wout  = squeeze(Wout_online_cell{m}(:,W_ind));
    % No transpose needed, online weights match stat
%     if m < 4
%         ww = reshape(Wout,128,128)';
%         Wout = ww(:);
%     end
    % Wout = [Wout;nan];
    for d = 1:num_days
        if ~isempty(bmi_weights_cell{m,d}) 
            [~,n_weights_mean] = calc_weights_from_mask_11042023(Wout,stat_cell{m,d});
            cur_b = bmi_weights_cell{m,d};
            cur_stds = z_stds{m,d}';
            
            % Check correlation with neural changes
            cur_means = all_overall_n_means{m,d};
            cur_diffs = cur_means(2,:)-cur_means(1,:);
            
            % Do I also do absolute changes. probably
            va_corrs(1,m,d) = corr(abs(cur_diffs'),abs(n_weights_mean));
            va_corrs(2,m,d) = corr(abs(cur_diffs'),abs(n_weights_mean).*cur_stds);
            va_corrs(3,m,d) = corr(abs(cur_diffs'),cur_stds);
            % va_corrs(m,d) = corr(cur_diffs',abs(n_weights_mean));
            
            for n = 1:nmodels
                % cur_n = squeeze(cur_b(n,1:end-1,yaw_ind));
                cur_n = -1.0*squeeze(cur_b(n,1:end-1,yaw_ind))';
                if m < 4
                    ww = reshape(cur_n,128,128)';
                    cur_n = ww(:);
                end
                cur_n = [cur_n;nan];
                [~,n_weights_mean_n] = calc_weights_from_mask_11042023(cur_n,stat_cell{m,d});
                % dot_prods(n,m,d) = sum((n_weights_mean_n.*n_weights_mean))/(norm(n_weights_mean_n)*norm(n_weights_mean));
                
                % Check correlation with neural changes
                yaw_corrs(1,n,m,d) = corr(abs(cur_diffs'),abs(n_weights_mean_n));
                yaw_corrs(2,n,m,d) = corr(abs(cur_diffs'),abs(n_weights_mean_n).*cur_stds);
                yaw_corrs(3,n,m,d) = corr(abs(cur_diffs'),cur_stds);
                % yaw_corrs(n,m,d) = corr(cur_diffs',abs(n_weights_mean_n));
            end
                
            
        end
    end
end