function [] = run_save_offline_2021_mice_bayes(mice, data_days, model_params, reps, training_trials, save_output)
% function for running offline analysis for many mice and days, and saving
% results (for plotting).
% Specifically for new cohort, DW113 ...

res_path = '/Users/ethan/Dropbox (Cambridge University)/Boston_PPC_Data_05_2021/'; 
% res_path = '/Users/ethansorrell/Dropbox (Cambridge University)/Boston_PPC_Data_05_2021/'; 
% train_test = 'test';
orig_model_params = model_params;
for m = 1:length(mice)
    data_d = data_days{m};
    for d = 1:length(data_d)
        
        model_params = orig_model_params;
        
        cur_res_path = res_path + mice(m) + '/' + data_d(d) + '/';
        xfull = importdata(cur_res_path + 'trimmed_virmen_' + mice(m) + '_' + data_d(d) + '_rearranged.mat');
        tbt_details = importdata(cur_res_path + mice(m) + '_' + data_d(d) + '_trial_by_trial_details.mat');
        t_summary = importdata(cur_res_path + mice(m) + '_' + data_d(d) + '_trial_types_summary.mat');
        if isfile(cur_res_path + mice(m) + '_' + data_d(d) + '_final.mat')
            zfull = importdata(cur_res_path + mice(m) + '_' + data_d(d) + '_final.mat');
        else
            z1 = importdata(cur_res_path + mice(m) + '_' + data_d(d) + '_final_1.mat');
            z2 = importdata(cur_res_path + mice(m) + '_' + data_d(d) + '_final_2.mat');
            zfull = [z1;z2];
        end
        
        zfull = permute_images(zfull);
        
        [ztrain, ztest, xtrain, xtest] = split_session_data(zfull, xfull, tbt_details, [1,2], training_trials, true, false);
        
        [ztrainfilt, train_mean, Rfixed] = preprocess_images_05_2021(ztrain,model_params,[],[],[],true);
        
        [ztestfilt, train_mean, Rfixed] = preprocess_images_05_2021(ztest,model_params,train_mean,Rfixed,[],false);
        
        model_params.dff = false;
        model_params.reg_images = false;
        model_params.spatial = false; 
        model_params.processed = true;
        
        for i = 1:reps
        
            % [results_struct] = DW_train_test_offline_func_05_2021(ztrainfilt, ztestfilt, xtrain, xtest, model_params, tbt_details, [1,2], 7, 4);
            [results_struct] = DW_train_test_offline_bayes(ztrainfilt, ztestfilt, xtrain, xtest, model_params, tbt_details,[1,2], 2, 3);
            
            if save_output
                % save('BMI_paper_figure_making/offline_results/res_struct_'+ mice(m) + '_' + data_d(d) + '_rep' + i + '_'+ training_trials + '.mat','results_struct');
                save('offline_results_bayes_alt_prior/res_struct_'+ mice(m) + '_' + data_d(d) + '_rep' + i + '_'+ training_trials + '.mat','results_struct');
            end
            
        end
        
    end
end
     