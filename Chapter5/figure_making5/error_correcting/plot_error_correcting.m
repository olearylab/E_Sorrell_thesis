%% Plot error correcting

addpath("../../../shared_functions")


virmen_cell = importdata("../../../shared_files/virmen_cell_5.mat");
tbt_cell = importdata("../../../shared_files/tbt_cell_5.mat");
virmen_train_cell = importdata("../../../shared_files/virmen_train_cell_5.mat");
tbt_train_cell = importdata("../../../shared_files/tbt_train_cell_5.mat");
centres = importdata("centres_50.mat"); % Change if different number of bins used

ex_md = [1,1];
offsets = [1.4953,1.4961;1.4953,1.4961;1.4953,1.4961;1.4804,1.4844;1.4804,1.4844];
linearise_x = true;
nbins = 50;

virmen_data = virmen_cell{ex_md};
tbt_details = tbt_cell{ex_md};

% plot_error_cartoon(virmen_data,tbt_details,nbins,offsets(1,:));

plot_specific_behaviour_examples_05072023(virmen_cell{1,2},tbt_cell{1,2},nbins,offsets(1,:),centres);
[mean_binned,std_binned] = calculate_mean_ball_va(virmen_train_cell{1},tbt_train_cell{1},nbins,offsets(1,:));
plot_specific_behaviour_examples_norm(virmen_cell{1,2},tbt_cell{1,2},nbins,offsets(1,:),mean_binned,std_binned,centres);

error_thresh = 1;
% [error_cell,x_cell,ro_all,sign_matching,h_boots] = run_many_check_error_correction_08052023(virmen_cell,tbt_cell,nbins);
[error_cell,x_cell,ro_all,sign_opposite,h_boots] = run_many_check_error_correction_norm(virmen_cell,tbt_cell,virmen_train_cell,tbt_train_cell,nbins,error_thresh);

% [error_cell,x_cell,ro_all,sign_matching,h_boots_dec] = run_many_check_BMI_correct_errors(virmen_cell,tbt_cell,nbins,true);
[error_cell,x_cell,ro_all,sign_opposite,h_boots_dec] = run_many_check_BMI_correct_errors_norm(virmen_cell,tbt_cell,virmen_train_cell,tbt_train_cell,nbins,error_thresh,true);


