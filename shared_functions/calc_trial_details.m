function [] = calc_trial_details(xfull,mouse,data_day)

% This function saves trial_by_trial_details and trial_types_summary

xfull(7,:) = wrapToPi(xfull(7,:));

% test_valid = xfull(8,:) == 0;
ITI = xfull(8,:);

% minor corrections
cleaned_valid = clean_valid_data(ITI);
% zdim = size(zdata,2);
trial_num_full = xfull(12,:);
% make sure ITI assigned to previous trial
trial_num_full = reassign_ITI_trial_num(cleaned_valid,trial_num_full);
trial_num_clean = trial_num_full(cleaned_valid);
num_trials = max(trial_num_clean);

trial_correct = xfull(9,:);
trial_correct_full = zeros(size(trial_correct));
trial_num = trial_num_full;
trial_type = xfull(1,:);
trial_lengths = zeros(num_trials,1);
trial_by_trial_correct = zeros(num_trials,1);
trial_by_trial_type = zeros(num_trials,1);

for j = 1:num_trials
    cur_trial_type = trial_type(trial_num == j);
    cur_trial_type = cur_trial_type(cleaned_valid(trial_num==j));
    cur_trial_xpos = xfull(5,trial_num==j);
    cur_trial_xpos = cur_trial_xpos(cleaned_valid(trial_num==j));
    trial_lengths(j) = sum(cleaned_valid(trial_num==j));
    trial_by_trial_type(j) = cur_trial_type(end);
    if sum(trial_correct(trial_num == j) ~= 0)
        trial_correct_full(trial_num == j) = 1;
        trial_by_trial_correct(j) = 1;
    elseif ismember(-1,ITI(trial_num == j))
        if (((trial_by_trial_type(j) == 1) || (trial_by_trial_type(j) == 3)) && cur_trial_xpos(end) < 0)
            trial_correct_full(trial_num == j) = 1;
            trial_by_trial_correct(j) = 1;
        elseif (((trial_by_trial_type(j) == 2) || (trial_by_trial_type(j) == 4)) && cur_trial_xpos(end) > 0)
            trial_correct_full(trial_num == j) = 1;
            trial_by_trial_correct(j) = 1;
        else
        trial_correct_full(trial_num == j) = -1;
        trial_by_trial_correct(j) = -1;
        end
    end
end

% Legacy for beahvioural statistics, might want to add back in though

% % Create summary statistics v2 file with structure:
% % Columns: mnlc, vnlc, mnlw, vnlw, mnlt, vnlt, mnrc, vnrc, mnrw, vnrw,
% % mnrt, vnrt, mblc, vblc, mblw, vblw, mblt, vblt, mbrc, vbrc, mbrw,
% % vbrw, mbrt, vbrt
% % rows: xpos, ypos, va, xvel, yvel, vavel, pitch, yaw, dpitch, dva,
% % triallength, pitch squared error
% virmen_summary_stats_v2 = zeros(12,24);
% xnums_ind = [5,6,7,2,3,4,13,15,16,17];
% % xnums_ind = [5,6,7,2,3,4,13,15];
% for j = 1:4
%     xfull_means = mean(xfull(:,test_valid&(trial_correct_full==1)&(trial_type==j)),2);
%     xfull_vars = var(xfull(:,test_valid&(trial_correct_full==1)&(trial_type==j)),[],2);
%     pitch_se = (xfull(13,test_valid&(trial_correct_full==1)&(trial_type==j)) - xfull(16,test_valid&(trial_correct_full==1)&(trial_type==j))).^2;
%     virmen_summary_stats_v2(1:10,6*(j-1)+1) = xfull_means(xnums_ind);
%     virmen_summary_stats_v2(11,6*(j-1)+1) = mean(trial_lengths((trial_by_trial_correct==1)&(trial_by_trial_type == j)));
%     virmen_summary_stats_v2(12,6*(j-1)+1) = mean(pitch_se);
%     virmen_summary_stats_v2(1:10,6*(j-1)+2) = xfull_vars(xnums_ind);
%     virmen_summary_stats_v2(11,6*(j-1)+2) = var(trial_lengths((trial_by_trial_correct==1)&(trial_by_trial_type == j)));
%     virmen_summary_stats_v2(12,6*(j-1)+2) = var(pitch_se);
% 
%     xfull_means = mean(xfull(:,test_valid&(trial_correct_full==-1)&(trial_type==j)),2);
%     xfull_vars = var(xfull(:,test_valid&(trial_correct_full==-1)&(trial_type==j)),[],2);
%     pitch_se = (xfull(13,test_valid&(trial_correct_full==-1)&(trial_type==j)) - xfull(16,test_valid&(trial_correct_full==-1)&(trial_type==j))).^2;
%     virmen_summary_stats_v2(1:10,6*(j-1)+3) = xfull_means(xnums_ind);
%     virmen_summary_stats_v2(11,6*(j-1)+3) = mean(trial_lengths((trial_by_trial_correct==-1)&(trial_by_trial_type == j)));
%     virmen_summary_stats_v2(12,6*(j-1)+3) = mean(pitch_se);
%     virmen_summary_stats_v2(1:10,6*(j-1)+4) = xfull_vars(xnums_ind);
%     virmen_summary_stats_v2(11,6*(j-1)+4) = var(trial_lengths((trial_by_trial_correct==-1)&(trial_by_trial_type == j)));
%     virmen_summary_stats_v2(12,6*(j-1)+4) = var(pitch_se);
% 
%     xfull_means = mean(xfull(:,test_valid&(trial_correct_full==0)&(trial_type==j)),2);
%     xfull_vars = var(xfull(:,test_valid&(trial_correct_full==0)&(trial_type==j)),[],2);
%     pitch_se = (xfull(13,test_valid&(trial_correct_full==0)&(trial_type==j)) - xfull(16,test_valid&(trial_correct_full==0)&(trial_type==j))).^2;
%     virmen_summary_stats_v2(1:10,6*(j-1)+5) = xfull_means(xnums_ind);
%     virmen_summary_stats_v2(11,6*(j-1)+5) = mean(trial_lengths((trial_by_trial_correct==0)&(trial_by_trial_type == j)));
%     virmen_summary_stats_v2(12,6*(j-1)+5) = mean(pitch_se);
%     virmen_summary_stats_v2(1:10,6*(j-1)+6) = xfull_vars(xnums_ind);
%     virmen_summary_stats_v2(11,6*(j-1)+6) = var(trial_lengths((trial_by_trial_correct==0)&(trial_by_trial_type == j)));
%     virmen_summary_stats_v2(12,6*(j-1)+6) = var(pitch_se);
% end

% Also save vector containing correct, wrong, and timeouts for each
% trial type
% format: nlc, nlw, nlt, nrc, nrw, nrt, blc, blw, blt, brc, brw, brt
trial_types_summary = zeros(12,1);
for j = 1:4
    trial_types_summary(3*(j-1)+1) = sum((trial_by_trial_type==j)&(trial_by_trial_correct==1));
    trial_types_summary(3*(j-1)+2) = sum((trial_by_trial_type==j)&(trial_by_trial_correct==-1));
    trial_types_summary(3*(j-1)+3) = sum((trial_by_trial_type==j)&(trial_by_trial_correct==0));
end

save(mouse+'/'+data_day+'/' +mouse+'_'+data_day+'_trial_types_summary.mat','trial_types_summary')

% Also save matrix with trial_by_trial_type, correct and vector of
% length equal to number of trials containing 
% number from 1-12 representing trial type + success corresponding to:
% nlc,nlw,nlt,nrc,nrw,nrt,blc,blw,blt,brc,brw,brt

trial_by_trial_details = zeros(3,length(trial_by_trial_type));

trial_by_trial_details(1,:) = trial_by_trial_type;
trial_by_trial_details(2,:) = trial_by_trial_correct;

trial_by_trial_details(3,trial_by_trial_correct == 1) = (trial_by_trial_type(trial_by_trial_correct == 1)-1)*3+1;
trial_by_trial_details(3,trial_by_trial_correct == -1) = (trial_by_trial_type(trial_by_trial_correct == -1)-1)*3+2;
trial_by_trial_details(3,trial_by_trial_correct == 0) = (trial_by_trial_type(trial_by_trial_correct == 0)-1)*3+3;

save(mouse+'/'+data_day+'/' +mouse+'_'+data_day+'_trial_by_trial_details.mat','trial_by_trial_details')

%correct_vec = [BMI_correct, num_trials_BMI, norm_correct, num_trials_norm, BMI_correct + norm_correct, num_trials];
%left_vs_right = [BMI_3_correct, BMI_3, BMI_4_correct, BMI_4, norm_1_correct, norm_1, norm_2_correct, norm_2];

% Legacy saved files
% correct_vec = [trial_types_summary(7)+trial_types_summary(10),sum(trial_types_summary(7:end)),trial_types_summary(1)+trial_types_summary(4),sum(trial_types_summary(1:6)),sum(trial_types_summary([1,4,7,10])),sum(trial_types_summary)];
% left_vs_right = [trial_types_summary(7),sum(trial_types_summary(7:9)),trial_types_summary(10),sum(trial_types_summary(10:end)),trial_types_summary(1),sum(trial_types_summary(1:3)),trial_types_summary(4),sum(trial_types_summary(4:6))];
% 
% save(cur_path+'correct_vec.mat','correct_vec')
% save(cur_path+'left_vs_right.mat','left_vs_right')

