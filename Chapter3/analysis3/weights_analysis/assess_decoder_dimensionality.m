function [full_all_res] = assess_decoder_dimensionality(results_struct,model_params,ztest,xtest,tbt_details,xnum)
% 13/10/2022

% assess decoder performance using increasing number of predictors (in
% descending order from largest weight).

% Start with just the bias term

% Might need to change, assumes all 8 in Wout
xnums = [6,3,13,15,5,2,7,4];
model_params.xnums = xnums(xnum);
t_types = [1,2];
if xnum == 7
    va_ind = 1;
else
    va_ind = [];
end

% Preprocess data
% initialise_params;
if ndims(ztest) == 2
    model_params.spatial = false;
    model_params.reg_images = false;
else
    model_params.reg_images = true;
    model_params.spatial = true;
end
create_reg = false;
disp("Pre-Processing")
[zdatafilt, train_mean, Rfixed] = preprocess_images_05_2021(ztest,model_params,results_struct.train_mean,results_struct.Rfixed,[],create_reg);
disp("Pre-Processing Complete")
model_params.spatial = false;
model_params.reg_images = false;
model_params.dff = false;

zdatafilt = [zdatafilt,ones(size(zdatafilt,1),1)];

Wout_full = results_struct.Wout(:,xnum);
results_struct.Wout = Wout_full;
results_struct.Wout(1:end-1,:) = 0;

full_all_res = zeros(3,size(Wout_full,1));

% [full_all_res(:,1)] = DW_test_only_07092022(zdatafilt,xtest,tbt_details,model_params,results_struct,[1,2], va_ind, []);
% predict results
xprediction = zdatafilt*results_struct.Wout;

test_ITI = xtest(8,:);
[test_valid] = clean_valid_data(test_ITI);
xfull = xtest;
xtest = xtest';
xtest = xtest(:,model_params.xnums);

results_struct.xtest = xtest;
results_struct.xprediction = xprediction;
results_struct.test_valid = test_valid;

[RMSE, R_square, r_p, test_valid_new, xtest_new, xprediction_new] = calc_errors_R_05_2021(results_struct,xfull,tbt_details,t_types,va_ind,[],[]);
full_all_res(:,1) = [RMSE;R_square;r_p];
% results_struct.all_res = all_res;

[sorted_w,w_inds] = sort(abs(Wout_full(1:end-1)),'descend');
for i = 1:(size(Wout_full,1)-1)
    results_struct.Wout(w_inds(i)) = Wout_full(w_inds(i));
    % [full_all_res(:,i+1)] = DW_test_only_07092022(zdatafilt,xtest,tbt_details,model_params,results_struct,[1,2], va_ind, []);
    xprediction = zdatafilt*results_struct.Wout;
    results_struct.xprediction = xprediction;
    [RMSE, R_square, r_p, test_valid_new, xtest_new, xprediction_new] = calc_errors_R_05_2021(results_struct,xfull,tbt_details,t_types,va_ind,[],[]);
    full_all_res(:,i+1) = [RMSE;R_square;r_p];
end
