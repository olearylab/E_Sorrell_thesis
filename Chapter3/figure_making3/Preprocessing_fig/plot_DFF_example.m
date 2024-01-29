function [] = plot_DFF_example(zdata,res_struct_raw,res_struct_dff)
% 25/05/2023

% Plotting example of raw neural data, and neural data post DF/F filtering.
% Also plot of the decoder output at the same points.

% Will need to be 4x2 - plotting early session and late session as difficult
% to show the whole transition.

figure

ax1 = subplot(2,2,1);
plot(movmean(squeeze(zdata(:,100,100)),1000))
% plot(squeeze(zdata(:,100,100)))

initialise_params;
model_params.reg_images = false;
model_params.spatial = false;
[ztrainfilt, train_mean, Rfixed] = preprocess_images_05_2021(squeeze(zdata(:,100,100)),model_params,[],[],[],false);

ax2 = subplot(2,2,3);
plot(movmean(ztrainfilt,1000))
% plot(ztrainfilt)

ax3 = subplot(2,2,2);
plot(res_struct_raw.xtest(:,1))
hold on
plot(res_struct_raw.xprediction(:,1))

ax4 = subplot(2,2,4);
plot(res_struct_dff.xtest(:,1))
hold on
plot(res_struct_dff.xprediction(:,1))

linkaxes([ax1,ax2,ax3,ax4],'x');