function [z_binned,error_mat,z_error_corrs,z_error_peak_corrs] = compare_neural_to_error(zdata,virmen_data,tbt_details,model_params,normalise_z,nbins,offsets,plot_res)
% 03/04/2023

% check whether there is a correlation between neural activity amplitude
% and error signal

% Bin by position first then check correlation.
% Don't 'need' to bin by position really though?

% error mat is trials x bins
[error_mat,x_binned] = check_error_correction_07032023(virmen_data,tbt_details,nbins,offsets,plot_res);

types_vec = [1,4,7,10];
max_types = max(tbt_details(1,:));
types_vec = types_vec(1:max_types);
%% remove invalid data
ITI = virmen_data(8,:);
cleaned_valid = clean_valid_data(ITI);

virmen_data = virmen_data(:,cleaned_valid);
%% Preprocess and bin data

initialise_params;
if ndims(zdata) == 2
    model_params.spatial = false;
    create_reg = false;
else
    create_reg = true;
end
model_params.reg_images = false;
disp("Pre-Processing")
[zdatafilt, train_mean, Rfixed] = preprocess_images_05_2021(zdata,model_params,[],[],cleaned_valid,create_reg);
disp("Pre-Processing Complete")   

%% Normalise (z-score) post filtering, before bootstrapping.
% Normalise filtered data with invalid data removed
if normalise_z 
    zdatafilt = zscore(zdatafilt);
end

xnum = 6;
linearise_x = true;
% z_binned is trials x nbins x zdim
[z_binned, centres] = bin_neural_data_general_any(virmen_data,zdatafilt,xnum,linearise_x,nbins);
zdim = size(z_binned,3);
% Not subsampling (for now at least)
sub_sample = false;
boot_samps = 100;
[bootstrap_means] = calc_bootsrapped_means_subsample(z_binned,tbt_details,boot_samps,types_vec,sub_sample);

z_binned_means = zeros(size(z_binned,2),size(z_binned,3),length(types_vec));
  
for i = 1:length(types_vec)
    z_binned_means(:,:,i) = mean(squeeze(bootstrap_means(i,:,:,:)),'omitnan');
end

[max_z_means,z_peaks_mean] = max(z_binned_means);

z_peaks_mean = squeeze(z_peaks_mean);
max_z_means = squeeze(max_z_means);
[~,z_large] = max(max_z_means(:,[3,4]),[],2);

z_error_corrs = nan.*ones(zdim,3);
z_error_peak_corrs = nan.*ones(zdim,3);
for n = 1:zdim
    cur_z = squeeze(z_binned(:,:,n));
    for i = 1:2    
        cur_trials = tbt_details(3,:)==types_vec(i+2);
        cur_zt = cur_z(cur_trials,:);
        cur_err = abs(error_mat(cur_trials,:));
        cur_corr = corrcoef(cur_zt(:),cur_err(:));
        z_error_corrs(n,i) = cur_corr(1,2);
        % peak
        cur_zpt = cur_zt(:,z_peaks_mean(n,i+2));
        cur_errp = cur_err(:,z_peaks_mean(n,i+2));
        cur_corr = corrcoef(cur_zpt(:),cur_errp(:));
        z_error_peak_corrs(n,i) = cur_corr(1,2);
    end
    cur_trials = ismember(tbt_details(3,:),types_vec([3,4]));
    cur_zt = cur_z(cur_trials,:);
    cur_err = abs(error_mat(cur_trials,:));
    cur_corr = corrcoef(cur_zt(:),cur_err(:));
    z_error_corrs(n,3) = cur_corr(1,2);   
    
    % peak
    cur_zpt = cur_zt(:,z_peaks_mean(n,z_large(n)+2));
    cur_errp = cur_err(:,z_peaks_mean(n,z_large(n)+2));
    cur_corr = corrcoef(cur_zpt(:),cur_errp(:));
    z_error_peak_corrs(n,3) = cur_corr(1,2);
end

if plot_res
    
    figure
    for i = 1:3
        scatter(i.*ones(zdim,1),z_error_corrs(:,i),'filled','k')
        hold on

    end
    yline(0,'--','LineWidth',2);
    xlim([0,4])
    
    figure
    for i = 1:3
        scatter(i.*ones(zdim,1),z_error_peak_corrs(:,i),'filled','k')
        hold on

    end
    yline(0,'--','LineWidth',2);
    xlim([0,4])
end

