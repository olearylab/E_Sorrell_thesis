function [z_binned_cell,error_cell,z_error_corrs_cell,z_error_peak_corrs_cell] = run_many_compare_neural_to_error(z_cell,virmen_cell,tbt_cell,model_params,normalise_z,nbins)
% 03/04/2023

% 
offsets = [1.4953,1.4961;1.4953,1.4961;1.4953,1.4961;1.4804,1.4844;1.4804,1.4844];
num_mice = size(z_cell,1);
num_days = size(z_cell,2);

mean_corrs = nan.*ones(num_mice,num_days);
std_corrs = nan.*ones(num_mice,num_days);

z_binned_cell = cell(num_mice,num_days);
error_cell = cell(num_mice,num_days);
z_error_corrs_cell = cell(num_mice,num_days);
z_error_peak_corrs_cell = cell(num_mice,num_days);

figure
for m = 1:num_mice
    for d = 1:num_days
        if ~isempty(z_cell{m,d})
            
            [z_binned,error_mat,z_error_corrs,z_error_peak_corrs] = compare_neural_to_error(z_cell{m,d},virmen_cell{m,d},tbt_cell{m,d},model_params,normalise_z,nbins,offsets(m,:),false);
            
            mean_corrs(m,d) = mean(z_error_corrs(:,3),'omitnan');
            std_corrs(m,d) = std(z_error_corrs(:,3),'omitnan');
            
            z_binned_cell{m,d} = z_binned;
            error_cell{m,d} = error_mat;
            z_error_corrs_cell{m,d} = z_error_corrs;
            z_error_peak_corrs_cell{m,d} = z_error_peak_corrs;
            
            for i = 1:3
                subplot(1,3,i)
                scatter(((m-1)*num_days+d)*ones(size(z_error_corrs,1),1),z_error_corrs(:,i),'k')
                % scatter(((m-1)*num_days+d)*ones(size(z_error_peak_corrs,1),1),z_error_peak_corrs(:,i),'k')
                hold on
            end
                
        end
    end
end

%% Plotting
for i = 1:3
    subplot(1,3,i)
    yline(0,'--','LineWidth',2);
end

figure
errorbar(mean_corrs(:),std_corrs(:))