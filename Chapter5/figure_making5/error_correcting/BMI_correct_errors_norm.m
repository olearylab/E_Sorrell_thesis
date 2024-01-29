function [error_mat,x_binned] = BMI_correct_errors_norm(virmen_data,tbt_details,mean_binned,std_binned,nbins,offsets)
% 09/09/2023

% Function for assessing how often the BMI output is error correcting. Just
% like checking if behaviour is error correcting.

% Start with just same bin/Sample, but maybe change to look at behaviour in
% previous bins/samples.

% Maybe smooth, maybe don't

% Work entirely binned for now (potentially change to samples).

% Look at whether it matches the ball yaw. Then maybe reduce to subset of
% error correcting ball yaw?

%%
% x_binned is trials x bins x xvec [6,7,15];
[error_mat,x_binned] = check_error_correction_norm(virmen_data,tbt_details,mean_binned,std_binned,nbins,offsets,false);

sign_matching = sign(squeeze(x_binned(:,:,3)))~=sign(error_mat);

types_vec = [1,4,7,10];

circum = 64;
V = 0.32;
alpha = -50/75*circum/V;
beta = 0.05*circum/V/2.5;
virmen_data(13,:) = alpha*(virmen_data(13,:)-offsets(1));
virmen_data(15,:) = -beta*(virmen_data(15,:)-offsets(2));

dt = 1/30;
% Differentiate view angle
vav = [virmen_data(17,2:end) - virmen_data(17,1:end-1),0]./dt;
virmen_data = [virmen_data;vav];
vav_ind = size(virmen_data,1);

% Low pass filter the velocities
% virmen_data(15,:) = lowpass(virmen_data(15,:),0.5,30);
% virmen_data(end,:) = lowpass(virmen_data(end,:),0.5,30);

linearise_x = true;
color_vec = [[0.4940 0.1840 0.5560];[0.4660 0.6740 0.1880]];
x_vec = [6,7,15,17,vav_ind];
virmen_data(7,:) = wrapToPi(virmen_data(7,:));
[x_binned,centres] = bin_kin_data(virmen_data,x_vec,linearise_x,nbins);
