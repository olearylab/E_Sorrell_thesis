%% Run offline neuron decoding

% Load in offline z_cell, virmen_cell, tbt_cell.
% Create the name vec of the mice+date
% set save_date to current date

initialise_params;

run_save_offline_neurons_15112022(z_cell, virmen_cell, tbt_cell, name_vec, model_params, 5, 80, true, save_date);