function [xpred, zreg] = Decoder_Online_Mod(zimage, Win, downsampled, afast, aslow, blurlarge, blursmall, dff, spatial, reg_images, train_mean, Rfixed)
%downsampling of image
if ~downsampled
    %[zvar] = downsample_image(zimage,4);
    zvar = imresize(zimage,1/4);
    zvar = zvar(:)';
else
    zvar = zimage(:)';
end

zvar = double(zvar);

% If registering images
if reg_images
    test_image = reshape(zvar',128,128);
    tformEstimate= imregcorr(test_image,train_mean,'translation');
    movingReg = imwarp(test_image,tformEstimate,'OutputView',Rfixed);
    zreg = movingReg;
    zvar = movingReg(:)';
else
    zreg = zimage;
end

%DF/F filtering
if dff
    [zvar] = dff_filt(zvar+150, afast, aslow);
end

%Spatial bandpass filtering
if spatial
    [zvar] = spatial_bandpass(zvar,blursmall,blurlarge);
end
    zvar = [zvar,1];

%Predict parameters and perform LMS (if update is true)
% This one for converting view angle to yaw
%[xpred]  = decoder_step(zvar, Win);
% This one for direct decoding
%[xpred]  = decoder_step_basic(zvar, Win);
% This one for low pass filtering view angle
[xpred] = decoder_step_lp(zvar,Win);

%[xpred, Wout]  = Predict_and_update(zvar, xvars, Win, lrate, valid, update);
%[xpred, Wout]  = Predict_and_update_ZA_LMS(zvar, xvars, Win, lrate, valid, update);