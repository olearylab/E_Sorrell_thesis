function [results_struct] = DW_train_test_offline_func_05_2021(ztrain, ztest, xtrain, xtest, model_params, tbt_details, t_types, va_ind, yaw_ind)


%%

% Instead decode sin and cos theta, then combine back to theta at the end
if model_params.va_2d
    sin_va = sin(xtrain(7,:));
    cos_va = cos(xtrain(7,:));
    xtrain = [xtrain; sin_va; cos_va; xtrain(7,:)];
    model_params.xnums = [model_params.xnums,size(xtrain,1)-2,size(xtrain,1)-1];
    
    % for including complex view angle
    % comp_va = exp(1i.*virmen_data(7,:));
    % virmen_data = [virmen_data; comp_va];
    % xnums = [xnums,size(virmen_data,1)];
    

end

tic
disp('training start')
[Wout, xpred, origztrainfilt, train_mean, Rfixed, train_valid, model_params] = DW_online_training_05_2021(ztrain, xtrain, model_params);
toc;
disp('training complete')
xtrainprediction = origztrainfilt*Wout;

xprediction = zeros(size(xtest,2),length(model_params.xnums));

% useful to not clear dff_filt between training and testing as we can't
% burn-in for testing.
%clear dff_filt;

ztest_reg = zeros(size(ztest));

tic;
disp('testing start')
for i = 1:size(ztest,1)
[xpred,zreg] = Decoder_Online_Mod_05_2021(squeeze(ztest(i,:,:)), Wout, model_params, train_mean,Rfixed);
xprediction(i,:) = xpred;
ztest_reg(i,:,:) = zreg;
end
toc;
disp('testing complete')

% If using sin and cos of view angle
if model_params.va_2d
    xprediction = [xprediction, atan2(xprediction(:,end-1),xprediction(:,end))];
    % model_params.xnums = model_params.xnums(1:end-2);
    % xprediction = xprediction(:,1:end-2);
    sin_va = sin(xtest(7,:));
    cos_va = cos(xtest(7,:));
    xtest = [xtest; sin_va; cos_va; xtest(7,:)];
    model_params.xnums = [model_params.xnums,size(xtest,1)];
    
end

test_ITI = xtest(8,:);
[test_valid] = clean_valid_data(test_ITI);
% test_valid = test_valid == 0;
% train_valid = xtrain(8,:);
% train_valid = train_valid == 0;

xtrain = xtrain';
xtrain = xtrain(:,model_params.xnums);
xfull = xtest;
xtest = xtest';
xtest = xtest(:,model_params.xnums);

results_struct.xtest = xtest;
results_struct.xtrain = xtrain;
results_struct.xprediction = xprediction;
results_struct.xtrainprediction = xtrainprediction;
results_struct.train_valid = train_valid;
results_struct.test_valid = test_valid;
results_struct.Wout = Wout;
results_struct.train_mean = train_mean;
results_struct.Rfixed = Rfixed;
results_struct.model_params = model_params;

% Come back to t_types, probably make an input. Also some of these others
% t_types = [1,2,3,4];
% va_ind = 7;
% yaw_ind = 4;
yaw_offset = 1.492; % Need to determine/ get off Dan.
[RMSE, R_square, r_p, test_valid_new, xtest_new, xprediction_new] = calc_errors_R_05_2021(results_struct,xfull,tbt_details,t_types,va_ind,yaw_ind,yaw_offset);
all_res = [RMSE;R_square;r_p];
results_struct.all_res = all_res;