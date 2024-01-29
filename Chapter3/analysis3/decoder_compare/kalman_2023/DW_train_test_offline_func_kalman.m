function [results_struct] = DW_train_test_offline_func_kalman(ztrain, ztest, xtrain, xtest, model_params, tbt_details, t_types, va_ind, yaw_ind)
% 26/05/2023
gain_weight = 1;
%%
xnums = model_params.xnums;
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

%%
tic
disp('training start')
[Wout, xpred, origztrainfilt, train_mean, Rfixed, train_valid, model_params] = DW_online_training_05_2021(ztrain, xtrain, model_params);
toc;
disp('training complete')
xtrainprediction = origztrainfilt*Wout;
xprediction = zeros(size(xtest,2),length(xnums));
xkfilt = zeros(size(xtest,2),length(xnums));

[C,A,R,Q,P] = kalman_trainer_2023(xtrain(xnums,:),origztrainfilt,Wout);

% clear dff_filt;

test_ITI = xtest(8,:);
[test_valid] = clean_valid_data(test_ITI);

tic;
disp('testing start')
xpost = xtest(xnums,1);
for i = 1:size(xtest,2)
% [xpred] = Decoder_Online_Mod(squeeze(ztest(i,:,:)), Wout, downsampled, afast, aslow, blurlarge, blursmall, dff, spatial);
[xpred,zreg] = Decoder_Online_Mod_05_2021(squeeze(ztest(i,:,:)), Wout, model_params, train_mean,Rfixed);
xprediction(i,:) = xpred;

    if test_valid(i) == 1
        if i>1
            if test_valid(i-1) ~= 1 
                xpost = xtest(xnums,i);
            end
        end
        [xpost, P] = kalman_filt_2023(C,A,R,Q,P,xpred,xpost,gain_weight);
        xkfilt(i,:) = xpost;
    else
        xkfilt(i,:) = nan;
        xkfilt(i,:) = nan;
    end
end
toc;
disp('testing complete')

xprediction = xkfilt;

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

% test_ITI = xtest(8,:);
% [test_valid] = clean_valid_data(test_ITI);
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