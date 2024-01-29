function [ss,centres] = prep_data_for_neural_bvn_shuffle(zdatafilt,virmen_data,tbt_details,nbins,keep_incorrect,kept_types,balance_types)
% 02/08/2023

% Prepare data for input into LSTM function.

types_vec = [1,4,7,10];
virmen_data(7,:) = wrapToPi(virmen_data(7,:));
% [x_binned,centres] = bin_kin_data(virmen_data,1:size(virmen_data,1),true,nbins);
%% Preprocess and bin data

% initialise_params;
% if ndims(zdata) == 2
%     model_params.spatial = false;
%     create_reg = false;
% else
%     create_reg = true;
% end
% model_params.reg_images = false;
% % disp("Pre-Processing")
% [zdatafilt, train_mean, Rfixed] = preprocess_images_05_2021(zdata,model_params,[],[],[],create_reg);
% % disp("Pre-Processing Complete")   

%% Normalise (z-score) post filtering, before bootstrapping.
% Normalise filtered data with invalid data removed
% if normalise_z 
%     zdatafilt = zscore(zdatafilt);
% end

xnum = 6;
linearise_x = true;
% z_binned is trials x nbins x zdim
[z_binned, centres] = bin_neural_data_general_any(virmen_data,zdatafilt,xnum,linearise_x,nbins);
zdim = size(z_binned,3);

if ~keep_incorrect
    z_binned = z_binned(tbt_details(2,:)==1,:,:);
else
    z_binned = z_binned(tbt_details(2,:)~=0,:,:);
end
z_binned_perm = permute(z_binned,[1,3,2]);
ss.binnedNeuralData2 = z_binned_perm;



trial_type = tbt_details(1,:);

if ~keep_incorrect
    trial_type = trial_type(tbt_details(2,:)==1);
    ss.correctVec = tbt_details(2,tbt_details(2,:)==1);
    tbt_details_cut = tbt_details(:,tbt_details(2,:)==1);
else
    trial_type = trial_type(tbt_details(2,:)~=0);
    ss.correctVec = tbt_details(2,tbt_details(2,:)~=0);
    tbt_details_cut = tbt_details(:,tbt_details(2,:)~=0);
end

ss.trialType = trial_type(ismember(trial_type,kept_types));

ss.correctVec = ss.correctVec(ismember(trial_type,kept_types));

ss.binnedNeuralData2 = ss.binnedNeuralData2(ismember(trial_type,kept_types),:,:);

tbt_details_cut = tbt_details_cut(:,ismember(trial_type,kept_types));

if balance_types
    if ~keep_incorrect
        % Should change to taking many subsamples, and averaging over these.
        sep_trial_nums = zeros(length(kept_types),1);
        for i = 1:length(kept_types)
            sep_trial_nums(i) = sum(tbt_details_cut(3,:) == types_vec(kept_types(i)));
        end
        kept_num = min(sep_trial_nums(sep_trial_nums~=0));

        all_kept = [];
        for i = 1:length(kept_types)
            cur_trials = find(tbt_details_cut(3,:)==types_vec(kept_types(i)));
            if ~isempty(cur_trials)
                cur_kept = datasample(cur_trials,kept_num,'Replace',false);
                all_kept = [all_kept,cur_kept];
            end
        end
    
        % Don't want to test on incorrect trials in this case?
        % inc_kept = find(tbt_details_cut(2,:)==-1);
        % all_kept = [all_kept,inc_kept];
    else
        
        % Should change to taking many subsamples, and averaging over these.
        sep_trial_nums = zeros(length(kept_types),1);
        for i = 1:length(kept_types)
            sep_trial_nums(i) = sum(ss.trialType == kept_types(i));
        end
        kept_num = min(sep_trial_nums(sep_trial_nums~=0));

        all_kept = [];
        for i = 1:length(kept_types)
            cur_trials = find(ss.trialType==kept_types(i));
            if ~isempty(cur_trials)
                cur_kept = datasample(cur_trials,kept_num,'Replace',false);
                all_kept = [all_kept,cur_kept];
            end
        end  
        
    end
    
    
    ss.binnedNeuralData2 = ss.binnedNeuralData2(all_kept,:,:);
    ss.trialType = ss.trialType(all_kept);
    ss.correctVec = ss.correctVec(all_kept);
    
end

% shuffle left and right trials separately
for i = 1:2
    cur_trials = find(ismember(ss.trialType,[1+(i-1),3+(i-1)]));
    shuffle_vec = randperm(length(cur_trials));

    ss.trialType(cur_trials) = ss.trialType(cur_trials(shuffle_vec));
end

ss.origTrialType = ss.trialType;

ss.trialType(ss.trialType==2) = 1;
ss.trialType(ss.trialType==3) = 2;
ss.trialType(ss.trialType==4) = 2;

% shuffle_vec = randperm(length(ss.trialType));
% ss.trialType = ss.trialType(shuffle_vec);
% ss.origTrialType = ss.origTrialType(shuffle_vec);

ss.correct_only = ~keep_incorrect;