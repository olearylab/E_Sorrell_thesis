function [Wout, xpred, ztrainfilt] = LMS_training_Mod_New(zimages, xvars, downsampled, afast, aslow, blurlarge, blursmall, lrate, dff, spatial, train_valid, loops, Win)
update = true;
train_length = size(xvars,1);
ztrainfilt = zeros(size(zimages,1),size(zimages,2)^2+1);
if length(Win)==1
    W = zeros(size(zimages,2)^2+1,size(xvars,2));
    W(end,:) = mean(xvars,1);
else
    W = Win;
end

% to burn in dff filter for training
for i = 1:size(xvars,1)
%downsampling of image data
if ~downsampled
    zimage = squeeze(zimages((end+1)-i,:,:));
    %[zvar] = downsample_image(zimage,4);
    zvar = imresize(zimage,1/4,'box');
    zvar = zvar(:)';
else
    zvar = zimages((end+1)-i,:);
end

% DF/F filtering
if dff
    [ztrainfilt(i,1:end-1)] = dff_filt(zvar, afast, aslow);
else
    ztrainfilt(i,1:end-1) = zvar;
end

end

for i = 1:size(xvars,1)
%downsampling of image data
if ~downsampled
    zimage = squeeze(zimages(i,:,:));
    %[zvar] = downsample_image(zimage,4);
    zvar = imresize(zimage,1/4,'box');
    zvar = zvar(:)';
else
    zvar = zimages(i,:);
end

% DF/F filtering
if dff
    [ztrainfilt(i,1:end-1)] = dff_filt(zvar, afast, aslow);
else
    ztrainfilt(i,1:end-1) = zvar;
end

% Spatial bandpass filter of data
if spatial
    [ztrainfilt(i,1:end-1)] = spatial_bandpass(ztrainfilt(i,1:end-1),blursmall,blurlarge);
end
ztrainfilt(i,end) = 1;
end

%origztrainfilt = ztrainfilt;

% Can add something here for removing first x trials due to dff burn in
% time

%

for i = 1:loops
    shuffle = randperm(train_length);
    %xvars = xvars(shuffle,:);
    %ztrainfilt = ztrainfilt(shuffle,:);
    %train_valid = train_valid(shuffle);

    for j = 1:train_length
        % LMS
        [xpred, Wout]  = Predict_and_update(ztrainfilt(shuffle(j),:), xvars(shuffle(j),:), W, lrate, train_valid(shuffle(j)), update);
        %[xpred, Wout]  = Predict_and_update_ZA_LMS(ztrainfilt(j,:), xvars(j,:), W, lrate, train_valid(j), update);
    end
end