%% Run training comparison

% Didn't make anything for running for all of them.

% Need to load in testing data, and filter testing data, along with virmen,
% tbt details. Also need to load in the initialisation weights (in the
% right transpose).

initialise_params;
loops_vec = [1,5,10,15];
hp_cut = 1;
model_params.xnums = 7;

[all_res_cell,hp_va_var_cell] = quick_check_script(ztrainfilt, ztestfilt, virmen_train, virmen_test,trial_by_trial_details,model_params, loops_vec, hp_cut, WinT);

% Then re-run for all the sessions. 
% only did the 1st testing session for m1-3