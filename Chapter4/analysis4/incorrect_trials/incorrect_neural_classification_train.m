function [Mdl_time,Mdl_bin,many_Mdl_bin,train_mean,Rfixed] = incorrect_neural_classification_train(xtrain, ztrain, tbt_train, linearise_x,nbins,c_type)
% 02/12/2021
% Supervised classification of neural data, to classify incorrect trials,
% trained on correct trials.

% Make input flexibly neurons or pixels. bin_neural_data preprocesses data
% assuming it is pixels at the moment.

% Dan suggested SVM or logistic regression

% Train various classifiers, then test on other days. Maybe training on
% first day testing on others isn't ideal (instead train on correct, or
% some of the correct trials from each day?). But will stick with this for
% now.

% Train classifiers for: all time points, all bins, separate bins

% Make input option of SVM or logistic.
% input needs to be 'svm' or 'logistic'.
% Still other classifiers to explore (and methods of computing classifier,
% and parameters to change).

%% Preprocess data
initialise_params;
if ndims(ztrain) == 2
    model_params.spatial = false;
    create_reg = false;
else
    create_reg = true;
end
model_params.reg_images = false;
disp("Pre-Processing")
[zdatafilt, train_mean, Rfixed] = preprocess_images_05_2021(ztrain,model_params,[],[],[],create_reg);
disp("Pre-Processing Complete")   


%% Bin neural data
% z_binned is trialsxbinsxneurons
% [z_binned_train, centres] = bin_neural_data(xtrain,ztrain',linearise_x,nbins);
[z_binned, centres] = bin_neural_data_general(xtrain,zdatafilt,linearise_x,nbins);

centres = centres*0.74;

%% All time points classifier

% remove invalid data
ITI = xtrain(8,:);
cleaned_valid = clean_valid_data(ITI);
% zdim = size(zdata,2);
trial_num_full = xtrain(12,:);
% make sure ITI assigned to previous trial
trial_num_full = reassign_ITI_trial_num(cleaned_valid,trial_num_full);
trial_num_clean = trial_num_full(cleaned_valid);
num_trials = max(trial_num_clean);

xtrain_clean = xtrain(:,cleaned_valid);
zdatafilt_clean = zdatafilt(cleaned_valid,:);

% label left turns as 1, right turns as 0. Only using correct trials for
% training so this is enough
train_labels = xtrain_clean(1,:) == 1;

% keep only correct trials
kept_trials = find(ismember(tbt_train(3,:),[1,4]));
zdatafilt_clean = zdatafilt_clean(ismember(trial_num_clean,kept_trials),:);
train_labels = train_labels(ismember(trial_num_clean,kept_trials))';

disp("Training Time Classifier")
Mdl_time = fitclinear(zdatafilt_clean,train_labels,'Learner',c_type);
disp("Training Complete")

%% all bins classifier
z_binned = z_binned(kept_trials,:,:);
z_binned_new = permute(z_binned,[3,2,1]);
z_binned_new = z_binned_new(:,:)';

bin_labels = tbt_train(1,:) == 1;
bin_labels = bin_labels(kept_trials);
bin_labels_new = repmat(bin_labels,[nbins,1]);
bin_labels_new = bin_labels_new(:);

disp("Training Bin Classifier")
Mdl_bin = fitclinear(z_binned_new,bin_labels_new,'Learner',c_type);
disp("Training Complete")

%% bin specific classifier
disp("Training Bin Specific Classifiers")
many_Mdl_bin = cell(nbins,1);
for i = 1:nbins
    many_Mdl_bin{i} = fitclinear(squeeze(z_binned(:,i,:)),bin_labels','Learner',c_type);
end
disp("Training Complete")
