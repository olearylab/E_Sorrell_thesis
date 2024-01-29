%% Run neural cue classifiers

% Don't have code for running all of them
% Load in training data then run 
linearise_x = true;
nbins = 50;
c_type = 'svm';
[Mdl_time,Mdl_bin,many_Mdl_bin,train_mean,Rfixed] = incorrect_neural_classification_train(xtrain, ztrain, tbt_train, linearise_x,nbins,c_type);

% Then load each testing set and run
[classifiers_results] = incorrect_neural_classification_test(ztest,xtest,linearise_x,nbins,Mdl_time,Mdl_bin,many_Mdl_bin,train_mean,Rfixed);