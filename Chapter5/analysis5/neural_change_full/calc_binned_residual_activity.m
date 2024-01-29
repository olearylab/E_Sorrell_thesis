function [z_binned] = calc_binned_residual_activity(zdata,virmen_data,tbt_details,model_params,normalise_z,nbins)
% 12/09/2023

%%
types_vec = [1,4,7,10];

%% remove invalid data
ITI = virmen_data(8,:);
cleaned_valid = clean_valid_data(ITI);

virmen_data = virmen_data(:,cleaned_valid);
%% Preprocess and bin data

if ndims(zdata) == 2
    model_params.spatial = false;
    create_reg = false;
else
    create_reg = true;
end
model_params.reg_images = false;
% disp("Pre-Processing")
[zdatafilt, train_mean, Rfixed] = preprocess_images_05_2021(zdata,model_params,[],[],cleaned_valid,create_reg);
% disp("Pre-Processing Complete")  

%% Normalise (z-score) post filtering, before bootstrapping.
% Normalise filtered data with invalid data removed
if normalise_z 
    [zdatafilt] = normalise_z_by_ball(zdatafilt,virmen_data,tbt_details);
end

xnum = 6;
linearise_x = true;
% z_binned is trials x nbins x zdim
[z_binned, centres] = bin_neural_data_general_any(virmen_data,zdatafilt,xnum,linearise_x,nbins);

% zdim = size(z_binned,3);
% Not subsampling (for now at least)
% Could just get mean as mean(z_binned(tbt_details(3,:) == types_vec(i),:,:));
sub_sample = false;
boot_samps = 100;
[bootstrap_means] = calc_bootsrapped_means_subsample(z_binned,tbt_details,boot_samps,types_vec,sub_sample);

z_binned_means = zeros(size(z_binned,2),size(z_binned,3),length(types_vec));
  
for i = 1:length(types_vec)
    z_binned_means(:,:,i) = mean(squeeze(bootstrap_means(i,:,:,:)),'omitnan');
end

mean_ind = [1,2,1,2];
for i = 1:length(types_vec)
    cur_trials = find(tbt_details(3,:) == types_vec(i));
    for n = 1:length(cur_trials)
        z_binned(cur_trials(n),:,:) = squeeze(z_binned(cur_trials(n),:,:)) - squeeze(z_binned_means(:,:,mean_ind(i)));
    end
end