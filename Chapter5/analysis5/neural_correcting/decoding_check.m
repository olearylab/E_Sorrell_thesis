function [all_res_mat,ball_weights,bmi_weights] = decoding_check(virmen_data,zdata,tbt_details,model_params,train_mean,Rfixed)

nModels = 5; % number of cv partitions
types_vec = [1,4,7,10];
kept_types = [1,2,3,4];
trial_num = virmen_data(12,:);

%% Preproccess neural data
[zfilt, train_mean, Rfixed] = preprocess_images_05_2021(zdata,model_params,train_mean,Rfixed,[],false);
model_params.dff = false;
model_params.spatial = false;
model_params.reg_images = false;

ball_weights = nan.*ones(nModels,size(zfilt,2)+1,length(model_params.xnums));
bmi_weights = nan.*ones(nModels,size(zfilt,2)+1,length(model_params.xnums));

%% Balance trial types (subsample)
% Should change to taking many subsamples, and averaging over these.
% Only choose trials that are correct and "good"
max_vas = nan.*ones(size(tbt_details,2),1);
for i = 1:size(tbt_details,2)
    % NB: Should change to only look at valid samples!!! Need to change
    % everywhere.
    max_vas(i) = max(abs(virmen_data(7,trial_num==i)));
end
tbt_details(3,max_vas>(2*pi)) = 13; % artificially set poor trials to 13.

sep_trial_nums = zeros(length(kept_types),1);
for i = 1:length(kept_types)
    sep_trial_nums(i) = sum(tbt_details(3,:) == types_vec(kept_types(i)));
end
kept_num = min(sep_trial_nums(sep_trial_nums~=0));

all_kept = [];
for i = 1:length(kept_types)
    cur_trials = find(tbt_details(3,:)==types_vec(kept_types(i)));
    if ~isempty(cur_trials)
        cur_kept = datasample(cur_trials,kept_num,'Replace',false);
        all_kept = [all_kept,cur_kept];
    end
end

virmen_data = virmen_data(:,ismember(trial_num,all_kept));
zfilt = zfilt(ismember(trial_num,all_kept),:);
sorted_kept = sort(all_kept);
orig_tbt_details = tbt_details;
tbt_details = tbt_details(:,sorted_kept);
trial_num = virmen_data(12,:);

%% Train on ball trials
ball_trials = [find(tbt_details(3,:) == types_vec(1)),find(tbt_details(3,:) == types_vec(2))];
% c = cvpartition(length(s.controlLabelSequence),'KFold',nModels);
% Choosing to repeat one split, rather than use 2, to balance left and
% right
c = cvpartition(length(ball_trials)/2,'KFold',nModels);
testSplitID = NaN(size(tbt_details,2));

ball_trial_nums = sorted_kept(ball_trials);

bmi_trials = [find(tbt_details(3,:) == types_vec(3)),find(tbt_details(3,:) == types_vec(4))];
bmi_trial_nums = sorted_kept(bmi_trials);

% for convenience save all in one.
% 1 = train ball/ test ball. 2 = train ball/ test bmi. 3 = train bmi/ test
% ball. 4 = train bmi/ test bmi
all_res_mat = nan.*ones(4,nModels,3,length(model_params.xnums));
%% Train network
for cvSplitNumber=1:nModels
    trainIdx = [c.training(cvSplitNumber);c.training(cvSplitNumber)];
    testIdx = [c.test(cvSplitNumber);c.test(cvSplitNumber)];

    XTrain = virmen_data(:,ismember(trial_num,ball_trial_nums(trainIdx)));
    ZTrain = zfilt(ismember(trial_num,ball_trial_nums(trainIdx)),:);

    [Wout, xpred, origztrainfilt, train_mean, Rfixed, train_valid, model_params] = DW_online_training_06_2023(ZTrain, XTrain, model_params);

    ball_weights(cvSplitNumber,:,:) = Wout;
    
    XTest = virmen_data(:,ismember(trial_num,ball_trial_nums(testIdx)));
    ZTest = zfilt(ismember(trial_num,ball_trial_nums(testIdx)),:);

    xprediction = zeros(size(XTest,2),length(model_params.xnums));

    for i = 1:size(ZTest,1)
    [xpred,zreg] = Decoder_Online_Mod_05_2021(squeeze(ZTest(i,:)), Wout, model_params, train_mean,Rfixed);
    xprediction(i,:) = xpred;
    end

    test_ITI = XTest(8,:);
    [test_valid] = clean_valid_data(test_ITI);
    % test_valid = test_valid == 0;
    % train_valid = xtrain(8,:);
    % train_valid = train_valid == 0;

    xtrain = XTrain';
    xtrain = XTrain(:,model_params.xnums);
    xfull = XTest;
    xtest = XTest';
    xtest = xtest(:,model_params.xnums);

    results_struct.xtest = xtest;
    results_struct.xtrain = xtrain;
    results_struct.xprediction = xprediction;
    % results_struct.xtrainprediction = xtrainprediction;
    % results_struct.train_valid = train_valid;
    results_struct.test_valid = test_valid;
    results_struct.Wout = Wout;
    results_struct.train_mean = train_mean;
    results_struct.Rfixed = Rfixed;
    results_struct.model_params = model_params;

    % Could create a fake tbt_details that says to test all trials as XTest is
    % only testing trials
    [RMSE, R_square, r_p, test_valid_new, xtest_new, xprediction_new] = calc_errors_R_05_2021(results_struct,XTest,orig_tbt_details,[1,2],7,[],[]);
    all_res_mat(1,cvSplitNumber,:,:) = [RMSE;R_square;r_p];

    bmi_XTest = virmen_data(:,ismember(trial_num,bmi_trial_nums));
    bmi_ZTest = zfilt(ismember(trial_num,bmi_trial_nums),:);

    xprediction = zeros(size(bmi_XTest,2),length(model_params.xnums));
    for i = 1:size(bmi_ZTest,1)
    [xpred,zreg] = Decoder_Online_Mod_05_2021(squeeze(bmi_ZTest(i,:)), Wout, model_params, train_mean,Rfixed);
    xprediction(i,:) = xpred;
    end

    test_ITI = bmi_XTest(8,:);
    [test_valid] = clean_valid_data(test_ITI);
    % test_valid = test_valid == 0;
    % train_valid = xtrain(8,:);
    % train_valid = train_valid == 0;

    xfull = bmi_XTest;
    xtest = bmi_XTest';
    xtest = xtest(:,model_params.xnums);

    results_struct.xtest = xtest;
    % results_struct.xtrain = xtrain;
    results_struct.xprediction = xprediction;
    % results_struct.xtrainprediction = xtrainprediction;
    % results_struct.train_valid = train_valid;
    results_struct.test_valid = test_valid;
    % results_struct.Wout = Wout;
    % results_struct.train_mean = train_mean;
    % results_struct.Rfixed = Rfixed;
    % results_struct.model_params = model_params;

    [RMSE, R_square, r_p, test_valid_new, xtest_new, xprediction_new] = calc_errors_R_05_2021(results_struct,bmi_XTest,orig_tbt_details,[3,4],7,[],[]);
    all_res_mat(2,cvSplitNumber,:,:) = [RMSE;R_square;r_p];

    % Not used
    testSplitID(testIdx)=cvSplitNumber;
end

%% Train on BMI trials
% bmi_trials = [find(tbt_details(3,:) == types_vec(3)),find(tbt_details(3,:) == types_vec(4))];
% bmi_trial_nums = sorted_kept(bmi_trials);
% c = cvpartition(length(s.controlLabelSequence),'KFold',nModels);
% Choosing to repeat one split, rather than use 2, to balance left and
% right
c = cvpartition(length(bmi_trials)/2,'KFold',nModels);
testSplitID = NaN(size(tbt_details,2));

% for convenience save all in one.
% 1 = train ball/ test ball. 2 = train ball/ test bmi. 3 = train bmi/ test
% ball. 4 = train bmi/ test bmi

%% Train network
for cvSplitNumber=1:nModels
    trainIdx = [c.training(cvSplitNumber);c.training(cvSplitNumber)];
    testIdx = [c.test(cvSplitNumber);c.test(cvSplitNumber)];

    XTrain = virmen_data(:,ismember(trial_num,bmi_trial_nums(trainIdx)));
    ZTrain = zfilt(ismember(trial_num,bmi_trial_nums(trainIdx)),:);

    [Wout, xpred, origztrainfilt, train_mean, Rfixed, train_valid, model_params] = DW_online_training_06_2023(ZTrain, XTrain, model_params);

    bmi_weights(cvSplitNumber,:,:) = Wout;
    
    XTest = virmen_data(:,ismember(trial_num,bmi_trial_nums(testIdx)));
    ZTest = zfilt(ismember(trial_num,bmi_trial_nums(testIdx)),:);

    xprediction = zeros(size(XTest,2),length(model_params.xnums));

    for i = 1:size(ZTest,1)
    [xpred,zreg] = Decoder_Online_Mod_05_2021(squeeze(ZTest(i,:)), Wout, model_params, train_mean,Rfixed);
    xprediction(i,:) = xpred;
    end

    test_ITI = XTest(8,:);
    [test_valid] = clean_valid_data(test_ITI);
    % test_valid = test_valid == 0;
    % train_valid = xtrain(8,:);
    % train_valid = train_valid == 0;

    xtrain = XTrain';
    xtrain = XTrain(:,model_params.xnums);
    xfull = XTest;
    xtest = XTest';
    xtest = xtest(:,model_params.xnums);

    results_struct.xtest = xtest;
    results_struct.xtrain = xtrain;
    results_struct.xprediction = xprediction;
    % results_struct.xtrainprediction = xtrainprediction;
    % results_struct.train_valid = train_valid;
    results_struct.test_valid = test_valid;
    results_struct.Wout = Wout;
    results_struct.train_mean = train_mean;
    results_struct.Rfixed = Rfixed;
    results_struct.model_params = model_params;

    % Could create a fake tbt_details that says to test all trials as XTest is
    % only testing trials
    [RMSE, R_square, r_p, test_valid_new, xtest_new, xprediction_new] = calc_errors_R_05_2021(results_struct,XTest,orig_tbt_details,[3,4],7,[],[]);
    all_res_mat(4,cvSplitNumber,:,:) = [RMSE;R_square;r_p];

    ball_XTest = virmen_data(:,ismember(trial_num,ball_trial_nums));
    ball_ZTest = zfilt(ismember(trial_num,ball_trial_nums),:);

    xprediction = zeros(size(ball_XTest,2),length(model_params.xnums));
    for i = 1:size(ball_ZTest,1)
    [xpred,zreg] = Decoder_Online_Mod_05_2021(squeeze(ball_ZTest(i,:)), Wout, model_params, train_mean,Rfixed);
    xprediction(i,:) = xpred;
    end

    test_ITI = ball_XTest(8,:);
    [test_valid] = clean_valid_data(test_ITI);
    % test_valid = test_valid == 0;
    % train_valid = xtrain(8,:);
    % train_valid = train_valid == 0;

    xfull = ball_XTest;
    xtest = ball_XTest';
    xtest = xtest(:,model_params.xnums);

    results_struct.xtest = xtest;
    % results_struct.xtrain = xtrain;
    results_struct.xprediction = xprediction;
    % results_struct.xtrainprediction = xtrainprediction;
    % results_struct.train_valid = train_valid;
    results_struct.test_valid = test_valid;
    % results_struct.Wout = Wout;
    % results_struct.train_mean = train_mean;
    % results_struct.Rfixed = Rfixed;
    % results_struct.model_params = model_params;

    [RMSE, R_square, r_p, test_valid_new, xtest_new, xprediction_new] = calc_errors_R_05_2021(results_struct,ball_XTest,orig_tbt_details,[1,2],7,[],[]);
    all_res_mat(3,cvSplitNumber,:,:) = [RMSE;R_square;r_p];

    % Not used
    testSplitID(testIdx)=cvSplitNumber;
end