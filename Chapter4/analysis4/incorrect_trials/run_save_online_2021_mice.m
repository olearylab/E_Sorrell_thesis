function [results_struct_cell] = run_save_online_2021_mice(mice, data_days, model_params,Wout_train_cell, train_mean_cell, Rfixed_cell, save_output)
% function for running offline analysis for many mice and days, and saving
% results (for plotting).
% Specifically for new cohort, DW113 ...

res_path = '/Users/ethan/Dropbox (Cambridge University)/Boston_PPC_Data_05_2021/'; 
% res_path = '/Users/ethansorrell/Dropbox (Cambridge University)/Boston_PPC_Data_05_2021/'; 
% train_test = 'test';
results_struct_cell = cell(5,4);
t_types = [1,2];
orig_model_params = model_params;
for m = 1:length(mice)
    data_d = data_days{m};
    res_struct.Wout = Wout_train_cell{m+3};
    res_struct.train_mean = train_mean_cell{m+3};
    res_struct.Rfixed = Rfixed_cell{m+3};
    for d = 1:length(data_d)
        
        model_params = orig_model_params;
        
        cur_res_path = res_path + mice(m) + '/' + data_d(d) + '/test/';
        xfull = importdata(cur_res_path + 'trimmed_virmen_' + mice(m) + '_' + data_d(d) + '_test_rearranged.mat');
        tbt_details = importdata(cur_res_path + mice(m) + '_' + data_d(d) + '_trial_by_trial_details.mat');
        t_summary = importdata(cur_res_path + mice(m) + '_' + data_d(d) + '_trial_types_summary.mat');
        zfull = importdata(cur_res_path + mice(m) + '_' + data_d(d) + '_final_test.mat');
        
        zfull = permute_images(zfull);
             
        [results_struct] = DW_test_only_20062023(zfull,xfull,tbt_details,model_params,res_struct,t_types, 7, 4);
        

        if save_output
            % save('BMI_paper_figure_making/offline_results/res_struct_'+ mice(m) + '_' + data_d(d) + '_rep' + i + '_'+ training_trials + '.mat','results_struct');
            % save('BMI_paper_figure_making/offline_results_raw/res_struct_'+ mice(m) + '_' + data_d(d) + '_rep' + i + '_'+ training_trials + '.mat','results_struct');
        end
        
        results_struct_cell{m+3,d} = results_struct;
        
    end
end
   