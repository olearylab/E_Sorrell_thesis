function [all_res_cell,hp_va_var_cell] = quick_check_script(ztrainfilt, ztestfilt, virmen_train, virmen_test,trial_by_trial_details,model_params, loops_vec, hp_cut, WinT)
    model_params.Win = WinT(:,2);
    model_params.remove_poor = false;
    [all_res_cell{1},hp_va_var_cell{1}] = check_overfitting_va_variance_effect_quick(ztrainfilt,ztestfilt,virmen_train,virmen_test,trial_by_trial_details,model_params, loops_vec, hp_cut, true);

    model_params.remove_poor = true;
    [all_res_cell{3},hp_va_var_cell{3}] = check_overfitting_va_variance_effect_quick(ztrainfilt,ztestfilt,virmen_train,virmen_test,trial_by_trial_details,model_params, loops_vec, hp_cut, true);

    model_params.Win = 0;
    [all_res_cell{4},hp_va_var_cell{4}] = check_overfitting_va_variance_effect_quick(ztrainfilt,ztestfilt,virmen_train,virmen_test,trial_by_trial_details,model_params, loops_vec, hp_cut, true);

    model_params.remove_poor = false;
    [all_res_cell{2},hp_va_var_cell{2}] = check_overfitting_va_variance_effect_quick(ztrainfilt,ztestfilt,virmen_train,virmen_test,trial_by_trial_details,model_params, loops_vec, hp_cut, true);

% 
% 1. = keep poor trials, initialise W
% 2. = keep poor trials, don't initialise W
% 3. = remove poor trials, initialise W
% 4. = remove poor trials, don't initialise W