%% Plot behave_change_ex

% load in usual virmen_cell and tbt_cell

nbins = 50;
linearise_x = true;
num_trials = 2;
offsets = [1.4953,1.4961;1.4953,1.4961;1.4953,1.4961;1.4804,1.4844;1.4804,1.4844];
plot_single_trials_ball(virmen_cell{1,1},tbt_cell{1,1},nbins,linearise_x,offsets(1,:),num_trials);
ex_md = [1,1];
plot_example_trajectories_06122022(virmen_cell, tbt_cell, ex_md, offsets(1,:), linearise_x, nbins);

% run the behavioural correlation code, or load in saved cells.
[h_boots] = plot_behav_corrs(behav_corr_cell,behav_std_corr_cell,behav_full_corr_cell);

% load in centres, and svm results
% plot_many_svm_bvn_classifier_results_18052023(accuracies_cell_bvn,centres)
[h_boots] = plot_many_svm_bvn_classifier_results_25072023(accuracies_cell_bvn,accuracies_cell_shuff,centres);