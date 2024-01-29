function [] = run_save_offline(mice, data_days, model_params, reps, training_trials, save_output)
% function for running offline analysis for many mice and days, and saving
% results (for plotting).
% Specifically for original cohort of DW81, 83 and 86
% Different function needed for new data.
res_path = '/Users/ethan/Dropbox (Cambridge University)/Boston_Data_09_2020/'; 
train_test = 'test';
orig_model_params = model_params;
for m = 1:length(mice)
    data_d = data_days{m};
    for d = 1:length(data_d)
        
        model_params = orig_model_params;
        
        cur_res_path = res_path + mice(m) + '/' + mice(m) + '_' + data_d(d) + '/' + train_test + '/';
        xfull = importdata(cur_res_path + 'trimmed_virmen_' + mice(m) + '_' + data_d(d) + '_' + train_test + '_revised.mat');
        tbt_details = importdata(cur_res_path + mice(m) + '_' + data_d(d) + '_trial_by_trial_details.mat');
        t_summary = importdata(cur_res_path + mice(m) + '_' + data_d(d) + '_trial_types_summary.mat');
        if isfile(cur_res_path + mice(m) + '_' + data_d(d) + '_final_test_revised.mat')
            zfull = importdata(cur_res_path + mice(m) + '_' + data_d(d) + '_final_test_revised.mat');
        else
            z1 = importdata(cur_res_path + mice(m) + '_' + data_d(d) + '_final_test_revised_1.mat');
            z2 = importdata(cur_res_path + mice(m) + '_' + data_d(d) + '_final_test_revised_2.mat');
            zfull = [z1;z2];
        end
        
        [ztrain, ztest, xtrain, xtest] = split_session_data(zfull, xfull, tbt_details, [1,2], training_trials, true, false);
        
        [ztrainfilt, train_mean, Rfixed] = preprocess_images_05_2021(ztrain,model_params,[],[],[],true);
        
        [ztestfilt, train_mean, Rfixed] = preprocess_images_05_2021(ztest,model_params,train_mean,Rfixed,[],false);
        
        model_params.dff = false;
        model_params.reg_images = false;
        model_params.spatial = false; 
        
        for i = 1:reps
        
            [results_struct] = DW_train_test_offline_func_05_2021(ztrainfilt, ztestfilt, xtrain, xtest, model_params, tbt_details, [1,2], 7, 4);
            
            if save_output
                % save('BMI_paper_figure_making/offline_results/res_struct_'+ mice(m) + '_' + data_d(d) + '_rep' + i + '_'+ training_trials + '.mat','results_struct');
                save('BMI_paper_figure_making/offline_results_raw/res_struct_'+ mice(m) + '_' + data_d(d) + '_rep' + i + '_'+ training_trials + '.mat','results_struct');
            end
            
        end
        
    end
end
        