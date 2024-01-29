function [z_binned,z_binned_res,x_binned,all_means] = calc_binned_residual_activity_from_samps(zdata,virmen_data,tbt_details,Wout,model_params,normalise_z,nbins,train_mean,Rfixed)
% 12/09/2023
% zdata can be neurons or pixels
%%
types_vec = [1,4,7,10];
virmen_data_full = virmen_data;
trial_num_full = virmen_data_full(12,:);

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
[zdatafilt_full, train_mean, Rfixed] = preprocess_images_05_2021(zdata,model_params,train_mean,Rfixed,[],create_reg);
zdatafilt = zdatafilt_full(cleaned_valid,:);
% disp("Pre-Processing Complete")  

%% Normalise (z-score) post filtering, before bootstrapping.
% Normalise filtered data with invalid data removed
% shouldn't need to normalise (as passing through decoder
if normalise_z 
    [zdatafilt] = normalise_z_by_ball(zdatafilt,virmen_data,tbt_details);
    [zdatafilt_full] = normalise_z_by_ball(zdatafilt_full,virmen_data_full,tbt_details);
end

xnum = 6;
linearise_x = true;
% z_binned is trials x nbins x zdim
[z_binned, centres] = bin_neural_data_general_any(virmen_data,zdatafilt,xnum,linearise_x,nbins);

virmen_data_full(6,:) = virmen_data_full(6,:) + abs(virmen_data_full(5,:));

% already cleaned above
virmen_data_clean = virmen_data_full(:,cleaned_valid);
[binned,edges] = discretize(virmen_data_clean(6,:),nbins);
[binned,edges] = discretize(virmen_data_full(6,:),edges);
% invalid data could be beyond final bin, so set to final bin
% This is just for making vav from decoder output sensible in the last
% sample of each trial.
binned(isnan(binned)) = nbins;
trial_num = virmen_data(12,:);

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

% Want to include ITI so when calculating vav don't introduce nans
z_res = NaN(size(zdatafilt_full));
mean_ind = [1,2,1,2];
for i = 1:length(types_vec)
    cur_trials = find(tbt_details(3,:) == types_vec(i));
    cur_trial_num = ismember(trial_num_full,cur_trials);
    for b = 1:nbins
        z_res((binned==b)&cur_trial_num,:) = zdatafilt_full((binned==b)&cur_trial_num,:) - squeeze(z_binned_means(b,:,mean_ind(i)));
    end
end

res_output = z_res*Wout;
on_output = zdatafilt_full*Wout;

dt = 1/30;
% Differentiate view angle
vav_res = [res_output(2:end,:) - res_output(1:end-1,:);0]./dt;
vav = [on_output(2:end,:) - on_output(1:end-1,:);0]./dt;

vav = vav(cleaned_valid);
vav_res = vav_res(cleaned_valid);

on_output = on_output(cleaned_valid);
res_output = res_output(cleaned_valid);

bmi_trials = find(ismember(tbt_details(3,:),types_vec([3,4])));
bmi_trial_num = ismember(trial_num,bmi_trials);
all_means = NaN(4,1);
all_means(1) = mean(abs(on_output(bmi_trial_num)));
all_means(2) = mean(abs(res_output(bmi_trial_num)));
all_means(3) = mean(abs(vav(bmi_trial_num)),'omitnan');
all_means(4) = mean(abs(vav_res(bmi_trial_num)),'omitnan');

[z_binned_res, centres] = bin_neural_data_general_any(virmen_data,z_res(cleaned_valid,:),xnum,linearise_x,nbins);

virmen_data = [virmen_data;on_output';res_output';vav';vav_res'];
on_ind = size(virmen_data,1)-3;
res_ind = size(virmen_data,1)-2;
vav_ind = size(virmen_data,1)-1;
vav_res_ind = size(virmen_data,1);

x_vec = [6,7,15,17,on_ind,res_ind,vav_ind,vav_res_ind];
virmen_data(7,:) = wrapToPi(virmen_data(7,:));
linearise_x = true; % not linearised in clean data yet
[x_binned,centres] = bin_kin_data(virmen_data,x_vec,linearise_x,nbins);