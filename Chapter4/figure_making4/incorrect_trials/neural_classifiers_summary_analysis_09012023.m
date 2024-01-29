function [h_boots] = neural_classifiers_summary_analysis_09012023(classifiers_res_cell,virmen_cell,tbt_cell,nbins)
% 06/12/2021

% Summarise results of svm/logistic regression classification.
% Plot classification accuracy for each trial type, and each mouse for each
% bin (or in time?)
% Potentially also plot output averaged across trials?
% Summarise results at cue offset, turn start, end of trial

% get and plot results for continuous time, bin, and bin specific classifiers

% Loop through several days. Store important results
num_days = size(virmen_cell,2);
num_mice = size(virmen_cell,1);

full_all_tbt = cell(num_mice,1);
% time results
full_all_res = cell(num_mice,1);
full_all_correct = cell(num_mice,1);
full_all_res_days = cell(num_mice,num_days);

% bin results
full_all_res_b = cell(num_mice,1);
full_all_correct_b = cell(num_mice,1);
full_all_res_days_b = cell(num_mice,num_days);

% bin specific results
full_all_res_bs = cell(num_mice,1);
full_all_correct_bs = cell(num_mice,1);
full_all_res_days_bs = cell(num_mice,num_days);

days_prop_bs = NaN(num_mice,num_days,nbins,4);

% for decoding turn direction
% l_trials = [1,5,7,11];
% r_trials = [2,4,8,10];
% for decoding cue direction
l_trials = [1,2,7,8];
r_trials = [4,5,10,11];
plot_types_t = [1,4;2,5;7,10;8,11];

for m = 1:num_mice
    cur_m_tbt = [];
    cur_m_res = [];
    cur_m_res_b = [];
    cur_m_res_bs = [];
    for d = 1:num_days
        if ~isempty(virmen_cell{m,d})
            %% Load in data
            classifiers_results = classifiers_res_cell{m,d};
            virmen_data = virmen_cell{m,d};
            tbt_details = tbt_cell{m,d};

            %% Neural check - time

            % Should I be checking match for each point then binning results,
            % or binning then checking. At the moment I am binning then
            % checking.
            % Possibly move in to separate function
            % x_binned is trials x bins x xdim
            % add classifier in time output to virmen data
            virmen_data = [virmen_data;classifiers_results{1}'];
            x_vec = size(virmen_data,1);
            linearise_x = true;
            [x_binned] = bin_kin_data(virmen_data,x_vec,linearise_x,nbins);

            cur_tbt = tbt_cell{m,d};
            cur_m_tbt = [cur_m_tbt,cur_tbt(3,:)];
            cur_m_res = [cur_m_res;squeeze(x_binned)];
            full_all_res_days{m,d} = squeeze(x_binned);

            %% Neural check - bins
            cur_m_res_b = [cur_m_res_b; reshape(classifiers_results{2},[nbins,length(classifiers_results{2})/nbins])'];
            full_all_res_days_b{m,d} = reshape(classifiers_results{2},[nbins,length(classifiers_results{2})/nbins])';
            %% Neural check - bin specific
            cur_m_res_bs = [cur_m_res_bs;classifiers_results{3}];
            cur_res_bs = classifiers_results{3};
            full_all_res_days_bs{m,d} = classifiers_results{3};
            
            %% For H_Boot
            % only works for completed trials, not timeouts
            cur_check = zeros(size(cur_tbt,2),1);
            cur_check(ismember(cur_tbt(3,:),l_trials)) = 1;
            cur_check(ismember(cur_tbt(3,:),r_trials)) = 0;

            cur_res_bs(cur_res_bs>=0.5) = 1;
            cur_res_bs(cur_res_bs<0.5) = 0;

            cur_correct_bs = cur_res_bs==repmat(cur_check,[1,nbins]);

            for i = 1:4
                days_prop_bs(m,d,:,i) = sum(cur_correct_bs(ismember(cur_tbt(3,:),plot_types_t(i,:)),:)==1,1)./sum(ismember(cur_tbt(3,:),plot_types_t(i,:)));
            end
        end
    end
    full_all_res{m} = cur_m_res;
    full_all_res_b{m} = cur_m_res_b;
    full_all_res_bs{m} = cur_m_res_bs;
    full_all_tbt{m} = cur_m_tbt;
    
    % only works for completed trials, not timeouts
    cur_m_check = zeros(length(cur_m_tbt),1);
    cur_m_check(ismember(cur_m_tbt,l_trials)) = 1;
    cur_m_check(ismember(cur_m_tbt,r_trials)) = 0;
    
    % ensure all values are 1 or 0
    cur_m_res(cur_m_res>=0.5) = 1;
    cur_m_res(cur_m_res<0.5) = 0;
    
    cur_m_res_b(cur_m_res_b>=0.5) = 1;
    cur_m_res_b(cur_m_res_b<0.5) = 0;
    
    cur_m_res_bs(cur_m_res_bs>=0.5) = 1;
    cur_m_res_bs(cur_m_res_bs<0.5) = 0;

    % repeat trial label for all bins. Check if matches decoder output.
    % result is 1 if match, 0 if not.
    full_all_correct{m} = cur_m_res==repmat(cur_m_check,[1,nbins]);
    full_all_correct_b{m} = cur_m_res_b==repmat(cur_m_check,[1,nbins]);
    full_all_correct_bs{m} = cur_m_res_bs==repmat(cur_m_check,[1,nbins]);
end

full_all_correct_3{1} = full_all_correct;
full_all_correct_3{2} = full_all_correct_b;
full_all_correct_3{3} = full_all_correct_bs;

%% Visualisation

% Get position of centre of each bin
virmen_ex = virmen_cell{1,1};
virmen_ex(6,:) = virmen_ex(6,:) + abs(virmen_ex(5,:));
[~,edges] = discretize(virmen_ex(6,virmen_ex(8,:)==0),nbins);

m_size = 50;
col_vec = ["b";"r";"b";"r"];
m_vec = ["o";"o";"^";"^"];
pos_scale = 0.74;
centres = pos_scale*(edges(2:end)+edges(1:end-1))/2;
edges = edges*pos_scale;
cue_end = 200*pos_scale;
turn_point = 300*pos_scale; % Could pick 305

plot_types_n = [1,2,4,5];
plot_types_b = [7,8,10,11];
titles_n = ["Normal Left Correct";"Normal Left Incorrect";"Normal Right Correct";"Normal Right Incorrect"];
titles_b = ["BMI Left Correct";"BMI Left Incorrect";"BMI Right Correct";"BMI Right Incorrect"];

% Loop through different classifiers, time, bin, bin specific
% for k = 1:3
%     cur_full_all_correct = full_all_correct_3{k};
%     final_results = zeros(m,nbins,8);
%     figure
%     for m = 1:num_mice
%         tbt_details = full_all_tbt{m};
%         cur_m_correct = cur_full_all_correct{m};
%         for j = 1:length(plot_types_n)
%             subplot(2,2,j)
%             hold on
%             xline(turn_point,'LineWidth',2,'Color','k');
%             xline(cue_end,'LineWidth',2,'Color','k');
%             yline(0,'LineWidth',2);
% 
%             cur_m_type_res = sum(cur_m_correct(tbt_details == plot_types_n(j),:)==1,1)./sum(tbt_details == plot_types_n(j));
%             final_results(m,:,j) = cur_m_type_res;
%             plot(centres,cur_m_type_res,'LineWidth',2)
% 
%             title(titles_n(j))
%             ylim([0,1])
%             if (j == 1) || (j == 3)
%                 ylabel("Classification Accuracy")
%                 yticks([0,0.5,1])
%             else
%                 yticks([])
%             end
%             if j > 2
%                 xlabel("Position (cm)")
%                 xticks([0,200])
%             else
%                 xticks([])
%             end
% 
% 
%         end
%     end
% 
%     figure
%     for m = 1:num_mice
%         tbt_details = full_all_tbt{m};
%         cur_m_correct = cur_full_all_correct{m};
%         for j = 1:length(plot_types_b)
%             subplot(2,2,j)
%             hold on
%             xline(turn_point,'LineWidth',2,'Color','k');
%             xline(cue_end,'LineWidth',2,'Color','k');
%             yline(0,'LineWidth',2);
% 
%             cur_m_type_res = sum(cur_m_correct(tbt_details == plot_types_b(j),:)==1,1)./sum(tbt_details == plot_types_b(j));
%             final_results(m,:,4+j) = cur_m_type_res;
%             plot(centres,cur_m_type_res,'LineWidth',2)
% 
%             title(titles_b(j))
%             ylim([0,1])
%             if (j == 1) || (j == 3)
%                 ylabel(["Classification"; "Accuracy"])
%                 yticks([0,0.5,1])
%             else
%                 yticks([])
%             end
%             if j > 2
%                 xlabel("Position (cm)")
%                 xticks([0,200])
%             else
%                 xticks([])
%             end
% 
% 
%         end
%     end
% 
% 
%     % Plot more condensed for cue_end, turn and end
%     figure
% 
%     % figure out which bins to keep
%     cue_bin = find(edges>=cue_end,1) - 1;
%     turn_bin = find(edges>=turn_point,1) - 1;
% 
%     titles = ["Cue End"; "Turn Start"; "End of Trial"];
%     bin_vec = [cue_bin,turn_bin,nbins];
%     for i = 1:3
%         subplot(1,3,i)
%         hold on
%         title(titles(i))
%         bar(mean(squeeze(final_results(:,bin_vec(i),:))),'w','LineWidth',2)
%         for j = 1:8
%             scatter(j*ones(num_mice,1),squeeze(final_results(:,bin_vec(i),j)),150,'k','filled')
%         end
%         if i == 1
%             ylabel("Classification Accuracy")
%             yticks([0,0.5,1])
%         else
%             yticks([])
%         end
%         ylim([0,1])
%         xlabel("Trial Type")
%         xticks(1:8)
%         xticklabels(["NLC";"NLI";"NRC";"NRI";"BLC";"BLI";"BRC";"BRI"])
%     end
% end

% Plot left and right collapsed
% Loop through different classifiers, time, bin, bin specific

plot_types_t = [1,4;2,5;7,10;8,11];
titles_t = ["Ball Correct";"Ball Incorrect";"BMI Correct";"BMI Incorrect"];

% for k = 1:3
% Only bin specific classifier
k = 3;
    cur_full_all_correct = full_all_correct_3{k};
    final_results = zeros(m,nbins,4);
    figure
    for m = 1:num_mice
        tbt_details = full_all_tbt{m};
        cur_m_correct = cur_full_all_correct{m};
        for j = 1:size(plot_types_t,1)
            subplot(2,2,j)
            hold on
            xline(turn_point,'--','LineWidth',2,'Color','k');
            xline(cue_end,'--','LineWidth',2,'Color','k');
            % yline(0,'LineWidth',2);

            cur_m_type_res = sum(cur_m_correct(ismember(tbt_details,plot_types_t(j,:)),:)==1,1)./sum(ismember(tbt_details,plot_types_t(j,:)));
            final_results(m,:,j) = cur_m_type_res;
            plot(centres,cur_m_type_res,'LineWidth',2,'Color',[0.5,0.5,0.5])

            title(titles_t(j))
            ylim([0,1])
            if (j == 1) || (j == 3)
                ylabel(["Classification"; "Accuracy"])
                yticks([0,0.5,1])
            else
                yticks([])
            end
            if j > 2
                xlabel("Position (cm)")
                xticks([0,200])
            else
                xticks([])
            end
            yline(0.5,'--','Linewidth',2);
            
        end
    end
    
    % Average lines
    mean_lines = squeeze(mean(final_results));
    
    for j = 1:size(plot_types_t,1)
        subplot(2,2,j)
        plot(centres,mean_lines(:,j),'LineWidth',3,'Color','k');
    end
    
    subplot(2,2,1)
    text(cue_end - 25, 0.2,["Cue";"End"],'FontSize',14,'HorizontalAlignment','center')
    text(turn_point - 25, 0.2,["Turn";"Start"],'FontSize',14,'HorizontalAlignment','center')
    % xticks([cue_end,turn_point])
    % xticklabels([{"Cue"; "End"};{"Turn"; "Start"}])

    % Plot more condensed for cue_end, turn and end
    %figure

    % figure out which bins to keep
    cue_bin = find(edges>=cue_end,1) - 1;
    turn_bin = find(edges>=turn_point,1) - 1;

%     titles = ["Cue End"; "Turn Start"; "End of Trial"];
%     bin_vec = [cue_bin,turn_bin,nbins];
%     for i = 1:3
%         subplot(1,3,i)
%         hold on
%         title(titles(i))
%         bar(mean(squeeze(final_results(:,bin_vec(i),:))),'w','LineWidth',2)
%         for j = 1:size(final_results,3)
%             scatter(j*ones(num_mice,1),squeeze(final_results(:,bin_vec(i),j)),150,'k','filled')
%         end
%         if i == 1
%             ylabel(["Classification"; "Accuracy"])
%             yticks([0,0.5,1])
%         else
%             yticks([])
%         end
%         ylim([0,1])
%         xlabel("Trial Type")
%         xticks(1:8)
%         xticklabels(["NC";"NI";"BC";"BI"])
%         yline(0.5,'LineWidth',2)
%     end
% end

%% Hierarchical Bootstraps
figure
% H bootstrap means and stds
rng(1);
boot_samps = 1000;
num_trials = 4;

% all_p_boot = NaN(nbins,1);
% for b = 1:nbins
%     [all_p_boot(b),all_centres(b,:),all_sems(b,:)] = run_H_boot_ets(squeeze(cur_mean(b,:,:)), squeeze(cur_mean_shuff(b,:,:)),false);
% end
orig_cur_means = days_prop_bs;
for m = 1:num_mice
    d_ind = 0;
    for d = 1:num_days
        if ~isnan(orig_cur_means(m,d,1,1)) 
            d_ind = d_ind+1;
            days_prop_bs(m,d_ind,:,:) = orig_cur_means(m,d,:,:);

        end
    end
    if d_ind<num_days
        days_prop_bs(m,d_ind+1:num_days,:,:) = nan;
    end
end

all_centres = NaN(nbins,4);
all_sems = NaN(nbins,4);
for b = 1:nbins
    for i = 1:4
        [bootstats] = get_bootstrapped_equalsamples(squeeze(days_prop_bs(:,:,b,i)),boot_samps,num_trials,'mean');
        all_centres(b,i) = mean(bootstats);
        all_sems(b,i) = std(bootstats);
    end
end
% h_boots.all_p_boot = all_p_boot;
h_boots.all_centres = all_centres;
h_boots.all_sems = all_sems;

% boostrap stds instead
boot_stds = all_sems;
lims_all = zeros(2,nbins,4);
lims_all(1,:,:) = all_centres - boot_stds;
lims_all(2,:,:) = all_centres + boot_stds;

for i = 1:size(plot_types_t,1)
    subplot(2,2,i)
    h = fill([centres,fliplr(centres)],[squeeze(lims_all(1,:,i)),fliplr(squeeze(lims_all(2,:,i)))],'k','EdgeColor','none');
    set(h,'facealpha',.3)
    hold on

    plot(centres,all_centres(:,i),'LineWidth',2,'Color','k')


% h = fill([centres,fliplr(centres)],[squeeze(lims_all(1,:,2)),fliplr(squeeze(lims_all(2,:,2)))],'k','EdgeColor','none');
% set(h,'facealpha',.3)
% hold on
% 
% plot(centres,all_centres(:,2),'--','LineWidth',2,'Color','k')


    yline(0.5,'--','Linewidth',2);
    xline(200*0.74,'--','Linewidth',2); % Cue end
    xline(300*0.74,'--','Linewidth',2); % Turn start
    % title(["Ball vs BMI Classifier";"Using Neural Data"])
    title(titles_t(i))
    ylim([0,1])
%     xlabel("Linearised Position (cm)")
% 
%     ylabel(["Classification Accuracy"])
    if (i == 1) || (i == 3)
        ylabel(["Classification"; "accuracy"])
        yticks([0,0.5,1])
    else
        yticks([])
    end
    if i > 2
        xlabel("Linearised position (cm)")
        xticks([0,200])
    else
        xticks([])
    end
    axis('square')
box off
end

subplot(2,2,1)
text(cue_end - 25, 0.2,["Cue";"end"],'FontSize',14,'HorizontalAlignment','center')
text(turn_point - 25, 0.2,["Turn";"start"],'FontSize',14,'HorizontalAlignment','center')