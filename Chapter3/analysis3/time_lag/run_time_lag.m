%% Run time lag decoding

% load in offline z_cell, virmen_cell, tbt_cell

lag_vec = [-30,-25,-20,-15,-10,-5,0,5,10,15,20,25,30];
[all_results] = compare_phase_lag_decoders_05062023(z_cell_CNN_offline, virmen_cell_offline, tbt_cell_offline, [1,2], true, lag_vec);