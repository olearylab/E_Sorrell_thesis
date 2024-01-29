%% Run assess decoder dimensionality

% The rest is all done in plotting code

% Load in the results_struct, ztest (pixels), xtest and tbt_details for the
% chosen session.

xnum = 1;
initialise_params;

[full_all_res] = assess_decoder_dimensionality(results_struct,model_params,ztest,xtest,tbt_details,xnum);