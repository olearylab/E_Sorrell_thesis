function [all_results] = compare_phase_lag_decoders_05062023(z_cell, virmen_cell, tbt_cell, t_types, neurons, lag_vec)
% 05/06/2022

% Editted from old function: compare_phase_lag_decoders_08_2022.

% function for assessing many decoder performances at different phase lags. 

% Use the same offline set of days

% input cell, one dimensional
% Need to edit if not using neurons

%%
initialise_params;

if neurons
    model_params.spatial = false;
    model_params.reg_images = false;
    model_params.lrate = 10^-2;
end

num_sesh = length(z_cell);

all_results = cell(num_sesh,1);

for i = 1:num_sesh
    zdata = z_cell{i};
    virmen_data = virmen_cell{i};
    tbt_details = tbt_cell{i};
    
    % results
    cur_results = cell(length(lag_vec),1);
        
    [ztrain, ztest, xtrain, xtest] = split_session_data(zdata, virmen_data, tbt_details, t_types, 80, true, false);
    
    % Ensure all the same training and testing behavioural samples used for
    % each lag.
    ITI_train = ~clean_valid_data(xtrain(8,:));
    ITI_test = ~clean_valid_data(xtest(8,:));

    neg_lag_train = [ITI_train(-lag_vec(1)+1:end),ones(1,-lag_vec(1))];
    pos_lag_train = [ones(1,lag_vec(end)),ITI_train(1:end-lag_vec(end))];

    ITI_train = ITI_train | pos_lag_train | neg_lag_train;

    neg_lag_test = [ITI_test(-lag_vec(1)+1:end),ones(1,-lag_vec(1))];
    pos_lag_test = [ones(1,lag_vec(end)),ITI_test(1:end-lag_vec(end))];

    ITI_test = ITI_test | pos_lag_test | neg_lag_test;

    xtrain(8,:) = ITI_train;
    xtest(8,:) = ITI_test;

    for l = 1:length(lag_vec)

    if lag_vec(l) < 0     

        ztrain_lag = ztrain(-lag_vec(l)+1:end,:,:);
        xtrain_lag = xtrain(:,1:end+lag_vec(l));
        ztest_lag = ztest(-lag_vec(l)+1:end,:,:);
        xtest_lag = xtest(:,1:end+lag_vec(l));
    else

        ztrain_lag = ztrain(1:end-lag_vec(l),:,:);
        xtrain_lag = xtrain(:,1+lag_vec(l):end);
        ztest_lag = ztest(1:end-lag_vec(l),:,:);
        xtest_lag = xtest(:,1+lag_vec(l):end);
    end

    % Train and test decoders

    [results_struct] = DW_train_test_offline_func_05_2021(ztrain_lag, ztest_lag, xtrain_lag, xtest_lag, model_params, tbt_details, [1,2], 7, 4);

    cur_results{l} = results_struct;

    end
    all_results{i} = cur_results;
end

%% plot results
% 
% plot_titles = ["ypos","yvel","pitch","yaw","xpos","xvel","VA","VAvel","VA from Yaw"];
% 
% for i = 1:3
%     figure(i)
%     for d = 1:length(data_days)
%         cur_plot = zeros(length(lag_vec),10);
%         for l = 1:length(lag_vec)
%             cur_res = all_results{d,l,2};
%             cur_plot(l,:) = cur_res(i,:);
%         end
%         for j = 1:9
%             subplot(3,3,j)
%             plot(lag_vec,cur_plot(:,j),'LineWidth',5)
%             title(plot_titles(j))
%             hold on
%         end
%     end
% end
