function [] = run_save_offline_neurons_15112022(z_cell, virmen_cell, tbt_cell, name_vec, model_params, reps, training_trials, save_output, save_date)
% 15/11/2022
% function for running offline analysis for many mice and days, and saving
% results (for plotting).

orig_model_params = model_params;

% As number of days per mouse can vary, just a 1d cell at the moment, name
% vec tells us what mice and session each one is. Maybe better solution
% (outer layer mice, inner layer days? that could work).

num_sessions = length(z_cell);


for d = 1:num_sessions

    model_params = orig_model_params;

    xfull = virmen_cell{d};
    tbt_details = tbt_cell{d};
    zfull = z_cell{d};

    [ztrain, ztest, xtrain, xtest] = split_session_data(zfull, xfull, tbt_details, [1,2], training_trials, true, false);

    [ztrainfilt, train_mean, Rfixed] = preprocess_images_05_2021(ztrain,model_params,[],[],[],false);

    [ztestfilt, train_mean, Rfixed] = preprocess_images_05_2021(ztest,model_params,train_mean,Rfixed,[],false);

    model_params.dff = false;
    model_params.reg_images = false;
    model_params.spatial = false; 

    for i = 1:reps

        [results_struct] = DW_train_test_offline_func_05_2021(ztrainfilt, ztestfilt, xtrain, xtest, model_params, tbt_details, [1,2], 7, 4);

        if save_output
            save('BMI_paper_figure_making/offline_neuron_results_CNN' + save_date + '/res_struct_'+ name_vec(d) + '_rep' + i + '_'+ training_trials + '.mat','results_struct');
        end

    end
end