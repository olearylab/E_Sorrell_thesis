function [Wout, xpred, origztrainfilt, train_mean, Rfixed] = DW_online_training(ztrain, virmen_data, xnums, lrate, loops, downsampled, dff, spatial, Win)
%zimages are downsampled
%binnedVirmen are the binned virmen data
%xnums are variables we are decoding

%trim virmen data
ztrain = ztrain(~isnan(virmen_data(1,:)),:,:);
virmen_data = virmen_data(:,~isnan(virmen_data(1,:)));
%wrap to pi view angle
virmen_data(7,:) = wrapToPi(virmen_data(7,:));

% For calculating mean training image for registration
% Should be one layer deeper so it's definitely on downsampled images.
if downsampled
    train_vec = double(ztrain(:,:));
    train_mean = mean(train_vec);
    train_mean = reshape(train_mean,128,128);
    Rfixed = imref2d(size(train_mean));
    clear train_vec
end

if dff
    ztrain = ztrain + 150; % add baseline to prevent negative values
end

xtrain = [];
for i = 1:length(xnums)
    xtrain = [xtrain, virmen_data(xnums(i),:)'];
end

train_length = size(xtrain,1);
trial_num = virmen_data(12,:);
invalid_sample = virmen_data(8,:);
train_invalid_sample = invalid_sample(1:train_length);
train_valid = train_invalid_sample == 0;

% Can comment out to keep incorrect trials for training
trial_correct = virmen_data(9,:);
train_trial_correct = trial_correct(1:train_length);
for i = 2:length(train_trial_correct)
    if trial_num(i) ~= trial_num(i-1)
        if train_trial_correct(i-1) == 0 
            train_valid(trial_num==trial_num(i-1)) = 0;
        end
    end
end

ztrain = double(ztrain);

fs = 30;

afast = (1-exp(-1/(0.5*fs)));
aslow = (1-exp(-1/(fs*45)));
%lrate = 10^-10;
blursmall = 0.6;
blurlarge = 5;

% remove samples where view angle is large. Changed to wrap to pi above
%train_valid(abs(virmen_data(7,:)')>2.5)=0;

%remove low and high values
% train_valid_old = train_valid;
% test_valid_old = test_valid;
% train_valid = train_valid&(xtrain(:,1)>5)'&(xtrain(:,1)<307)';
% test_valid = test_valid&(xtest(:,1)>5)'&(xtest(:,1)<307)';

% Ensure functions are reset
clear dff_filt
clear Predict_and_update

%% Training
%tic
%disp('training start')
[Wout, xpred, origztrainfilt] = LMS_training_Mod_New(ztrain, xtrain, downsampled, afast, aslow, blurlarge, blursmall, lrate, dff, spatial, train_valid, loops, Win);
%toc;
%disp('training complete')
