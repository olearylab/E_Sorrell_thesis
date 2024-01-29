%% Run offline position decoding

% this is what I used to get the results used for plotting. But I could
% change the plotting code to use results from more standard offline
% decoding. But not sure where that code is? - run_save_online and
% _2021_mice. using saved weights, and train means etc. Should check
% plotting with these results, and possibly replace. Should be simple
% matter of changing ind from 3 to 1. and loading in the other results
% struct (results_struct_cell_pixels) - just checked, results are
% identical.

% need to load training and testing data, and tbt_details
intialise_params;

[train_results_struct] = train_test_direction_decoder(ztrain, ztest, xtrain, xtest, model_params, tbt_details, [1,2]);

% Run after training so don't have to train separately each time.
[results_struct] = test_only_direction_decoder(ztest, xtest, model_params, tbt_details, [1,2], train_results_struct);