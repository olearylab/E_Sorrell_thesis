function [results_struct] = train_test_direction_decoder(ztrain, ztest, xtrain, xtest, model_params, tbt_details, t_types)
% 17/11/2021
%%
% for including complex view angle
% comp_va = exp(1i.*virmen_data(7,:));
% virmen_data = [virmen_data; comp_va];
% xnums = [xnums,size(virmen_data,1)];

% Set up to decode cue and view angle
% Arbitrarily set 2 as the direction (need to keep 1 for the training
% balancing)
% Set left as positive and right as negative, since this is the case in the
% normal data
% Also decode position (as a disengagement check)
% Maybe should just decode all the usual variables...
model_params.xnums = [2,7,6];
xtrain(model_params.xnums(1),xtrain(1,:)==1) = 1;
xtrain(model_params.xnums(1),xtrain(1,:)==2) = -1;

xtest(model_params.xnums(1),ismember(xtest(1,:),[1,3])) = 1;
xtest(model_params.xnums(1),ismember(xtest(1,:),[2,4])) = -1;

tic
disp('training start')
[Wout, xpred, origztrainfilt, train_mean, Rfixed, train_valid, model_params] = DW_online_training_05_2021(ztrain, xtrain, model_params);
toc;
disp('training complete')
xtrainprediction = origztrainfilt*Wout;

xprediction = zeros(size(xtest,1),length(model_params.xnums));

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
va_ind = 2;
yaw_ind = [];
yaw_offset = [];
% yaw_offset = 1.492; % Need to determine/ get off Dan.
[RMSE, R_square, r_p, test_valid_new, xtest_new, xprediction_new] = calc_errors_R_05_2021(results_struct,xfull,tbt_details,t_types,va_ind,yaw_ind,yaw_offset);
all_res = [RMSE;R_square;r_p];
results_struct.all_res = all_res;