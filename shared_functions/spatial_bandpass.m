function [bp_data] = spatial_bandpass(data, blursmall, blurlarge)
    zmat = reshape(data,128,128);
    zmat = imgaussfilt(zmat,blursmall)-imgaussfilt(zmat,blurlarge);
    bp_data = zmat(:)';