function [RMSE_all,R_square_all] = calc_many_online_decoder_errors(virmen_cell,tbt_cell)
% 31/05/2022

% Calculate control trial decoder errors.
% input cells are one d.

%% 
RMSE_all = zeros(length(virmen_cell),2);
R_square_all = zeros(length(virmen_cell),2);
for n = 1:length(virmen_cell)
    
    [RMSE, R_square, r_p, test_valid_new, xtest_new, xprediction_new] = calc_errors_R_online_results(virmen_cell{n},tbt_cell{n},[1,2],2,[],[],[13,7]);
    
    RMSE_all(n,:) = RMSE;
    R_square_all(n,:) = R_square;
    
end
    