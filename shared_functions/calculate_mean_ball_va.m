function [mean_binned,std_binned] = calculate_mean_ball_va(virmen_data,tbt_details,nbins,offsets)
% 09/09/2023

% calculate binned mean and std for view angle on training ball trials

types_vec = [1,4,7,10];
num_types = 2;

circum = 64;
V = 0.32;
alpha = -50/75*circum/V;
beta = 0.05*circum/V/2.5;
virmen_data(13,:) = alpha*(virmen_data(13,:)-offsets(1));
virmen_data(15,:) = -beta*(virmen_data(15,:)-offsets(2));

linearise_x = true;
% color_vec = [[0.4940 0.1840 0.5560];[0.4660 0.6740 0.1880]];
x_vec = [6,7,15,17]; % Changed from [6,7,15];
virmen_data(7,:) = wrapToPi(virmen_data(7,:));
[x_binned,centres] = bin_kin_data(virmen_data,x_vec,linearise_x,nbins);

% Old method, was calculating mean view angle from testing ball trials
mean_binned = nan.*ones(size(x_binned,2),size(x_binned,3),num_types);
std_binned = nan.*ones(size(x_binned,2),size(x_binned,3),num_types);
for i = 1:num_types
    mean_binned(:,:,i) = squeeze(mean(x_binned(tbt_details(3,:)==types_vec(i),:,:),1,'omitnan'));
    std_binned(:,:,i) = squeeze(std(x_binned(tbt_details(3,:)==types_vec(i),:,:),0,1,'omitnan'));
end