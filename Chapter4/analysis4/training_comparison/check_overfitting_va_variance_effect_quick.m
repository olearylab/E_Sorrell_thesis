function [all_res_mat,hp_va_var_mat] = check_overfitting_va_variance_effect_quick(ztrain,ztest,xtrain,xtest,tbt_details,model_params, loops_vec, hp_cut, is_filtered)
% 24/05/2022

% Function for checking the effect of overfitting on high-passed view angle
% variance
% Just for one mouse and day

% retrain decoder several times with different numbers of loops through
% training data
% Plot R^2 and hp view angle variance after each loop/ 5 loops
model_params.xnums = 7;
fs = 30;
control_trials = find(ismember(tbt_details(3,:),[1,4]));
%% Preprocess data once
if ~is_filtered
    disp("Pre-processing Training Data")
    model_params.reg_images = false;
    [ztrainfilt, train_mean, Rfixed] = preprocess_images_05_2021(ztrain,model_params,[],[],[],true);
    model_params.reg_images = true;
    disp("Pre-processing Testing Data")
    [ztestfilt, train_mean, Rfixed] = preprocess_images_05_2021(ztest,model_params,train_mean,Rfixed,[],false);
else 
    ztrainfilt = ztrain;
    ztestfilt = ztest;
end

model_params.dff = false;
model_params.reg_images = false;
model_params.spatial = false;

all_res_mat = zeros(3,length(loops_vec));

hp_va_var_mat = zeros(1,length(loops_vec));

% remove invalid data
ITI = xtest(8,:);
cleaned_valid = clean_valid_data(ITI);
xtest_clean = xtest(:,cleaned_valid);
trial_num = xtest_clean(12,:);

loops_vec(2:end) = loops_vec(2:end) - loops_vec(1:end-1);

for i = 1:length(loops_vec)
    model_params.loops = loops_vec(i);
    [results_struct] = DW_train_test_offline_func_05_2021(ztrainfilt, ztestfilt, xtrain, xtest, model_params, tbt_details, [1,2], 1, []);
    
    all_res_mat(:,i) = results_struct.all_res;
    
    all_processed = zeros(size(results_struct.xprediction,1),1);
    all_vars = zeros(size(tbt_details,2),1);
    
    % Low pass filter predicted view angle as with online data.
    xprediction_lp = online_lp_copy(results_struct.xprediction);
    xprediction_lp = xprediction_lp(cleaned_valid);
    
    for n = 1:size(tbt_details,2)
        cur_trial = xprediction_lp(trial_num==n);
        all_processed(trial_num==n) = highpass(cur_trial,hp_cut,fs);
        all_vars(n) = var(all_processed(trial_num==n));
    end
    
    hp_va_var_mat(i) = mean(all_vars(control_trials));
    model_params.Win = results_struct.Wout;
end

% %% Plotting
% 
% figure
% subplot(1,2,1)
% plot(loops_vec,hp_va_var_mat,'LineWidth',2)
% subplot(1,2,2)
% plot(loops_vec,all_res_mat(1,:),'LineWidth',2)