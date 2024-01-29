function [] = plot_many_compare_diff_training_12062023(full_all_res_cell,full_hp_va_var_cell)
% 12/06/2023

% plot results of checking effects of different things on training.

% mouse 1 2 and 3

% 1. = keep poor trials, initialise W
% 2. = keep poor trials, don't initialise W
% 3. = remove poor trials, initialise W
% 4. = remove poor trials, don't initialise W

% assume loops_vec is [1,5,10,15] - need to check this.

marker_size = 10;

num_mice = size(full_all_res_cell,1);

final_plots_res = zeros(num_mice,8);
final_plots_var = zeros(num_mice,8);

final_bars_res = zeros(num_mice,4,2);
final_bars_var = zeros(num_mice,4,2);

for i = 1:num_mice
    for j = 1:4
        cur_res = full_all_res_cell{i,j};
        cur_var = full_hp_va_var_cell{i,j};
        
        final_plots_res(i,2*(j-1)+1:2*j) = cur_res(1,[2,4]);
        final_plots_var(i,2*(j-1)+1:2*j) = cur_var([2,4]);
        
        final_bars_res(i,j,:) = cur_res(1,[2,4]);
        final_bars_var(i,j,:) = cur_var([2,4]);
        
    end
end

%% Plot

figure
subplot(1,2,1)
plot(final_plots_res','-o','LineWidth',2,'Color','k')

subplot(1,2,2)
plot(final_plots_var','-o','LineWidth',2,'Color','k')

%% Plot paired bars of 5/15 loops. gray and white bars?
x_shift = 0.15;
xx = [1,2,3,4];
x5 = xx-x_shift;
x15 = xx+x_shift;

figure
X = [1,2,3,4];
% X = reordercats(X,{'Small','Medium','Large','Extra Large'});
Y = squeeze(mean(final_bars_res));
bar(X,Y,'hist')
hold on
for m =1:num_mice
    plot(x5,squeeze(final_bars_res(m,:,1)),'o','Color','k','MarkerFaceColor','k','MarkerSize',marker_size)
    plot(x15,squeeze(final_bars_res(m,:,2)),'o','Color','k','MarkerFaceColor','k','MarkerSize',marker_size)
end

xlabel("Training Protocol")
ylabel(["Ball Trial";"View Angle RMSE (rad)"]) 
    

figure
X = [1,2,3,4];
% X = reordercats(X,{'Small','Medium','Large','Extra Large'});
Y = squeeze(mean(final_bars_var));
bar(X,Y,'hist')
hold on
for m =1:num_mice
    plot(x5,squeeze(final_bars_var(m,:,1)),'o','Color','k','MarkerFaceColor','k','MarkerSize',marker_size)
    plot(x15,squeeze(final_bars_var(m,:,2)),'o','Color','k','MarkerFaceColor','k','MarkerSize',marker_size)
end

xlabel("Training Protocol")
ylabel(["Ball Trial 1Hz HighPass";"View Angle Variance"]) 

%% Hierarchical bootstrap
