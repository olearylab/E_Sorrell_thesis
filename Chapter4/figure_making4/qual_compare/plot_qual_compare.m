%% plot qualitative comparisons

virmen_cell = importdata("../../../shared_files/virmen_cell_5.mat");
tbt_cell = importdata("../../../shared_files/tbt_cell_5.mat");

virmen_data = virmen_cell{1,1};
tbt_details = tbt_cell{1,1};
num_trials = 3;

plot_example_position_trajectories_thesis(virmen_data,tbt_details,num_trials);

nbins = 50;
linearise_x = true;
boot_samps = 100;
sub_sample = false;
ex_offset = 1.4804;

virmen_data_2 = virmen_cell{4,1};
tbt_details_2 = tbt_cell{4,1};

plot_pva_trajectory_comparisons_0612022(virmen_data_2,tbt_details_2,nbins,linearise_x,boot_samps,sub_sample,ex_offset);

[h_boots] = trial_length_comparisons_plot_thesis(virmen_cell,tbt_cell);