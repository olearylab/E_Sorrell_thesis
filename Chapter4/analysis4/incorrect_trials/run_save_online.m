function [results_struct_cell] = run_save_online(mice, data_days, model_params, Wout_train_cell, train_mean_cell, Rfixed_cell, save_output)
% 28/06/2023
% function for running decoding on online days for many mice and days, and saving
% results (for plotting).
% Specifically for original cohort of DW81, 83 and 86
% Different function needed for new data.
res_path = '/Users/ethan/Dropbox (Cambridge University)/Boston_Data_09_2020/'; 
% res_path = '/Users/ethansorrell/Dropbox (Cambridge University)/Boston_Data_09_2020/'; 
train_test = 'test';
orig_model_params = model_params;
results_struct_cell = cell(5,4);
t_types = [1,2];
for m = 1:length(mice)
    data_d = data_days{m};
    res_struct.Wout = Wout_train_cell{m};
    res_struct.train_mean = train_mean_cell{m};
    res_struct.Rfixed = Rfixed_cell{m};
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
        
        
        [results_struct] = DW_test_only_20062023(zfull,xfull,tbt_details,model_params,res_struct,t_types, 7, 4);
        

        if save_output
            % save('BMI_paper_figure_making/offline_results/res_struct_'+ mice(m) + '_' + data_d(d) + '_rep' + i + '_'+ training_trials + '.mat','results_struct');
            % save('BMI_paper_figure_making/offline_results_raw/res_struct_'+ mice(m) + '_' + data_d(d) + '_rep' + i + '_'+ training_trials + '.mat','results_struct');
        end
        
        results_struct_cell{m,d} = results_struct;
        
    end
end
      