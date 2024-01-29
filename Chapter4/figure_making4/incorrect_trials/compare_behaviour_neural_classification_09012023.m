function [conf_mat, full_compare, full_compare_correct] = compare_behaviour_neural_classification_09012023(classifiers_res_cell,results_struct_cell,virmen_cell,tbt_cell,nbins,start_dist,offset_vec)
% 07/12/21

% Function for comparing trials classification of incorrect trials using
% behaviour and neural activity

% Plot using some kind of Venn diagram

% combine all mice together? Get total number of incorrect trials
% classified correctly by only one method, and correctly by both methods.
% Then also get total number of incorrect trials classified incorrectly by
% one or the other or both.

% use start of turn bin for neural data. Use time classifier. Maybe change
% so this can be determined by input.

num_mice = size(virmen_cell,1);
num_days = size(virmen_cell,2);
pos_scale = 0.74;

virmen_ex = virmen_cell{1,1};
virmen_ex(6,:) = virmen_ex(6,:) + abs(virmen_ex(5,:));
[~,edges] = discretize(virmen_ex(6,virmen_ex(8,:)==0),nbins);

centres = pos_scale*(edges(2:end)+edges(1:end-1))/2;
edges = edges*pos_scale;
cue_end = 200*pos_scale;
turn_point = 300*pos_scale; % Could pick 305

turn_bin = find(edges>=turn_point,1) - 1;

%% Get relevant results
ends_match_combine_days = cell(num_mice,1);
full_m_tbt = cell(num_mice,1);
full_all_res = cell(num_mice,1);
full_all_correct = cell(num_mice,1);

% for decoding turn direction
% l_trials = [1,5,7,11];
% r_trials = [2,4,8,10];
% for decoding cue direction
l_trials = [1,2,7,8];
r_trials = [4,5,10,11];

full_compare = zeros(num_mice,4);
% also for correct trials
full_compare_correct = zeros(num_mice,4);
matching_trials = cell(num_mice,1);
non_matching_trials = cell(num_mice,1);

for m = 1:num_mice
    cur_m_tbt = [];
    cur_m_ends = [];
    cur_m_res = [];
    for d = 1:num_days
        if ~isempty(virmen_cell{m,d})
            [trials_mat,trial_ends,ends_match,trial_corrs] = end_yaw_integration(virmen_cell{m,d}, tbt_cell{m,d}, start_dist, offset_vec(m), false);

            cur_tbt = tbt_cell{m,d};
            cur_m_tbt = [cur_m_tbt,cur_tbt(3,:)];
            cur_m_ends = [cur_m_ends;ends_match];
            
            classifiers_results = classifiers_res_cell{m,d};
            virmen_data = virmen_cell{m,d};
            tbt_details = tbt_cell{m,d};
            
            virmen_data = [virmen_data;classifiers_results{1}'];
            x_vec = size(virmen_data,1);
            linearise_x = true;
            [x_binned] = bin_kin_data(virmen_data,x_vec,linearise_x,nbins);

            cur_m_res = [cur_m_res;squeeze(x_binned)];
        end
    end
    
    % convert ends match into 1 if match 0 if not
    % reverse as looking at if end matches cue
    % cur_m_ends(cur_m_ends==1) = 1;
    % cur_m_ends(cur_m_ends==-1) = 0;
    cur_m_ends(cur_m_ends==1) = 0;
    cur_m_ends(cur_m_ends==-1) = 1;
    
    full_m_tbt{m} = cur_m_tbt;
    ends_match_combine_days{m} = cur_m_ends;
    full_all_res{m} = cur_m_res;
    
    % only works for completed trials, not timeouts
    cur_m_check = zeros(length(cur_m_tbt),1);
    cur_m_check(ismember(cur_m_tbt,l_trials)) = 1;
    cur_m_check(ismember(cur_m_tbt,r_trials)) = 0;
    
    % ensure all values are 1 or 0
    cur_m_res(cur_m_res>=0.5) = 1;
    cur_m_res(cur_m_res<0.5) = 0;

    % repeat trial label for all bins. Check if matches decoder output.
    % result is 1 if match, 0 if not.
    full_all_correct{m} = cur_m_res==repmat(cur_m_check,[1,nbins]);
    
    % number behave only, neural only, both. Correct then incorrect classification.
    %cur_compare = zeros(2,3);
    cur_trials = ismember(cur_m_tbt,[8,11]);
    cur_m_correct = full_all_correct{m};
    full_compare(m,1) = sum((cur_m_ends(cur_trials) == 1) & (cur_m_correct(cur_trials,turn_bin) == 0));
    full_compare(m,2) = sum((cur_m_ends(cur_trials) == 0) & (cur_m_correct(cur_trials,turn_bin) == 1));
    full_compare(m,3) = sum((cur_m_ends(cur_trials) == 1) & (cur_m_correct(cur_trials,turn_bin) == 1));
    full_compare(m,4) = sum((cur_m_ends(cur_trials) == 0) & (cur_m_correct(cur_trials,turn_bin) == 0));
    
    matching_trials{m} = cur_trials(((cur_m_ends(cur_trials) == 1) & (cur_m_correct(cur_trials,turn_bin) == 1)) | ((cur_m_ends(cur_trials) == 0) & (cur_m_correct(cur_trials,turn_bin) == 0)));
    non_matching_trials{m} = cur_trials(((cur_m_ends(cur_trials) == 0) & (cur_m_correct(cur_trials,turn_bin) == 1)) | ((cur_m_ends(cur_trials) == 1) & (cur_m_correct(cur_trials,turn_bin) == 0)));
    
    % For correct trials
    cur_trials = ismember(cur_m_tbt,[7,10]);
    full_compare_correct(m,1) = sum((cur_m_ends(cur_trials) == 1) & (cur_m_correct(cur_trials,turn_bin) == 0));
    full_compare_correct(m,2) = sum((cur_m_ends(cur_trials) == 0) & (cur_m_correct(cur_trials,turn_bin) == 1));
    full_compare_correct(m,3) = sum((cur_m_ends(cur_trials) == 1) & (cur_m_correct(cur_trials,turn_bin) == 1));
    full_compare_correct(m,4) = sum((cur_m_ends(cur_trials) == 0) & (cur_m_correct(cur_trials,turn_bin) == 0));
end

%% Visualisation
% Use full_all_correct and ends_match_combine_days to get number in each
% category, then plot venn diagram, or venn style plot. Confusion matrix
% style
% figure
% 
conf_mat = [sum(full_compare(:,4)),sum(full_compare(:,2));sum(full_compare(:,1)),sum(full_compare(:,3))];
% 
% confusionchart(conf_mat,["Turn";"Cue"]);
% % confusionchart(conf_mat);
% title("Incorrect Trial Classifier Comparison")
% xlabel("Neural Classification")
% ylabel("Behavioural Classification")


% Plot showing proportion of incorrect trials with matching classification,
% as well as distributions of position decoding for matched vs unmatched
% trials.
% Actually just box plot of distributions.

%% TODO: Double check below works as expected

% get trial by trial accuracy
all_trials_cell = cell(num_mice,num_days);
all_trials_m_cell = cell(num_mice,1);
for m = 1:num_mice
cur_trials_m = [];
    for d = 1:num_days
        if ~isempty(virmen_cell{m,d})
            [all_trials_res,trials_summary_means,trials_summary_stds] = trial_by_trial_accuracy(results_struct_cell{m,d}, virmen_cell{m,d}, tbt_cell{m,d});
            all_trials_cell{m,d} = all_trials_res;
            cur_trials_m = [cur_trials_m;all_trials_res];
        end
    end
    all_trials_m_cell{m} = cur_trials_m;
end

ind = 3;
full_all_match = [];
full_all_no_match = [];
for m = 1:num_mice
    cur_trials = all_trials_m_cell{m};
    full_all_match = [full_all_match;squeeze(cur_trials(matching_trials{m},1,ind))];
    full_all_no_match = [full_all_no_match;squeeze(cur_trials(non_matching_trials{m},1,ind))];
end

combined_res = [full_all_match;full_all_no_match];
group_vec = [ones(length(full_all_match),1);2*ones(length(full_all_no_match),1)];

% Boxplot
% figure
% hold on
% boxplot(combined_res*pos_scale,group_vec)
% title("Position Decoding Accuracy")
% ylabel("RMSE (cm)")
% ylim([0,100*pos_scale])
% xlabel("Trial Type")
% box off