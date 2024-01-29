function [classifiers_results] = incorrect_neural_classification_test(ztest,xtest,linearise_x,nbins,Mdl_time,Mdl_bin,many_Mdl_bin,train_mean,Rfixed)
% 03/12/2021
% Supervised classification of neural data, to classify incorrect trials,
% trained on correct trials.

% Make input flexibly neurons or pixels. bin_neural_data preprocesses data
% assuming it is pixels at the moment.

% Dan suggested SVM or logistic regression

% For testing once we have trained models

% Test classifiers for: all time points, all bins, separate bins

% test all time points, remove invalid etc during plotting.

%% Preprocess data
initialise_params;
if ndims(ztest) == 2
    model_params.spatial = false;
    model_params.reg_images = false;
end
create_reg = false;
disp("Pre-Processing")
[zdatafilt, train_mean, Rfixed] = preprocess_images_05_2021(ztest,model_params,train_mean,Rfixed,[],create_reg);
disp("Pre-Processing Complete")   


%% Bin neural data
% z_binned is trialsxbinsxneurons
% [z_binned_train, centres] = bin_neural_data(xtrain,ztrain',linearise_x,nbins);
[z_binned, centres] = bin_neural_data_general(xtest,zdatafilt,linearise_x,nbins);

centres = centres*0.74;

%% All time points classifier

disp("Testing Time Classifier")
time_test = predict(Mdl_time,zdatafilt);
disp("Testing Complete")

%% all bins classifier
z_binned_new = permute(z_binned,[3,2,1]);
z_binned_new = z_binned_new(:,:)';

disp("Testing Bin Classifier")
bin_test = predict(Mdl_bin,z_binned_new);
disp("Testing Complete")

%% bin specific classifier
disp("Testing Bin Specific Classifiers")
bin_specific_test = zeros(size(z_binned,1),nbins);

for i = 1:nbins
    bin_specific_test(:,i) = predict(many_Mdl_bin{i},squeeze(z_binned(:,i,:)));
end
disp("Testing Complete")

classifiers_results{1} = time_test;
classifiers_results{2} = bin_test;
classifiers_results{3} = bin_specific_test;