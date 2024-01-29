function [results_struct] = DW_train_test_offline_bayes(ztrain, ztest, xtrain, xtest, model_params, tbt_details, t_types, va_ind, yaw_ind)
% 26/05/2023

%%
% Just worry about these 4, save time.
xnums = [6,7,13,15];
model_params.xnums = xnums;
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

nbins = 100;

if ~model_params.processed
    if ndims(ztrain) == 2
        model_params.spatial = false;
        create_reg = false;
    else
        create_reg = true;
    end
    orig_params = model_params;
    model_params.reg_images = false;
    disp("Pre-Processing")
    [ztrainfilt, train_mean, Rfixed] = preprocess_images_05_2021(ztrain,model_params,[],[],[],create_reg);
    disp("Pre-Processing Complete") 

    model_params.reg_images = orig_params.reg_images;
    disp("Pre-Processing")
    [ztestfilt, train_mean, Rfixed] = preprocess_images_05_2021(ztest,model_params,train_mean,Rfixed,[],false);
    disp("Pre-Processing Complete")  
else
    ztrainfilt = ztrain;
    ztestfilt = ztest;
end

%training
disp('Training start')
tic
%[meanyt,varpix,prior,Xedges,Yedges] = bayes_train(ztrain,xtrain,downsampled,dff,spatial,nbins);
% [meany,meanth,varpix,yprior,thprior,Yedges,Thedges] = bayes_separate_train(ztrain,xtrain,downsampled,dff,spatial,nbins);
[all_means,varpix,all_priors,all_edges] = bayes_separate_train(ztrainfilt,xtrain,nbins,xnums);
toc
disp('Training end')
%testing
xpred = zeros(size(ztest,1),length(xnums));

log_post_all = zeros(size(ztest,1),nbins,length(xnums));
log_like_all = zeros(size(ztest,1),nbins,length(xnums));
disp('Testing start')
tic
xprev = xtest(xnums,1);

test_ITI = xtest(8,:);
[test_valid] = clean_valid_data(test_ITI);

for i = 1:size(ztest,1)
    if test_valid(i) == 1
        if i>1
            if test_valid(i-1) ~= 1 
                xprev = xtest(xnums,i);
            end
        end
        %[log_post(i,:,:), ypred(i), thpred(i)] = bayes_test(squeeze(ztest(i,:,:)),prior,Xedges,Yedges,meanyt,varpix,downsampled,dff,spatial,nbins);
        % [log_like_all(i,:,:), log_post_all(i,:,:), xpred(i,:)] = bayes_separate_test(squeeze(ztest(i,:,:)),yprev,thprev,yprior,thprior,Yedges,Thedges,meany,meanth,varpix,downsampled,dff,spatial,nbins);
        [log_like_all(i,:,:), log_post_all(i,:,:), xpred(i,:)] = bayes_separate_test(ztestfilt(i,:),xprev,all_priors,all_edges,all_means,varpix,nbins); 
        xprev = xpred(i,:);
    else
        xpred(i,:) = nan;
    end
end
disp('Testing end')
toc
xprediction = xpred;

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
% results_struct.xtrainprediction = xtrainprediction;
% results_struct.train_valid = train_valid;
results_struct.test_valid = test_valid;
% results_struct.Wout = Wout;
% results_struct.train_mean = train_mean;
% results_struct.Rfixed = Rfixed;
results_struct.model_params = model_params;

% Come back to t_types, probably make an input. Also some of these others
% t_types = [1,2,3,4];
% va_ind = 7;
% yaw_ind = 4;
yaw_offset = 1.492; % Need to determine/ get off Dan.
[RMSE, R_square, r_p, test_valid_new, xtest_new, xprediction_new] = calc_errors_R_05_2021(results_struct,xfull,tbt_details,t_types,va_ind,yaw_ind,yaw_offset);
all_res = [RMSE;R_square;r_p];
results_struct.all_res = all_res;