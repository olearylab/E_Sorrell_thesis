function [] = plot_example_single_trajectories_16072023(virmen_cell, tbt_cell, ex_md, offsets, linearise_x, nbins, n_trials)
% 06/12/2022

% Plot example trajectories for pitch and yaw to compare between bmi and
% control trials.

% Maybe somehow combine with LSTM figures.
x_vec = [13,15];
%% 
% 2 x 2
% Best to collapse into 1 x 2 if possible
% Do I want actual time? Or binned?
types_vec = [1,4,7,10];

virmen_data = virmen_cell{ex_md(1),ex_md(2)};
tbt_details = tbt_cell{ex_md(1),ex_md(2)};

ITI = virmen_data(8,:);
cleaned_valid = clean_valid_data(ITI);
virmen_data = virmen_data(:,cleaned_valid);

% convert into velocity
circum = 64;
V = 0.32;
alpha = -50/75*circum/V;
beta = 0.05*circum/V/2.5;
virmen_data(13,:) = alpha*(virmen_data(13,:)-offsets(1));
virmen_data(15,:) = -beta*(virmen_data(15,:)-offsets(2));

% Bin data: trials x nbins x xdim
[x_binned,centres] = bin_kin_data(virmen_data,x_vec,linearise_x,nbins);

%% Plot averages and error bars - alternative
% Plot either mean or median line with shaded in standard deviation, or
% confidence interval, or something similar

% Pick better colour scheme?
plot_types = [1,7;4,10];
% Do I want lighter indivudal trials? Or lighter stds or something similar?
% Do I need to have separate colours in this case?
% plot_colour = [0 0.4470 0.7410];
plot_colours = [[0.4940 0.1840 0.5560];[0.4660 0.6740 0.1880]];
figure
pos_scale = .74;
centres = centres*pos_scale;

start_ind = 1;

% figure
for i = 1:2
    for j = 1:2 
        subplot(2,2,(i-1)*2+j)
        for t = 1:2
            cur_trials = find(ismember(tbt_details(3,:),types_vec(j+(t-1)*2)));
            for n = 1:n_trials
           
                plot(centres,squeeze(x_binned(cur_trials(start_ind+n),:,i)),'Color',plot_colours(t,:),'LineWidth',2)

                hold on
            end
        end
        box off
    end
end

ax11 = subplot(2,2,1);
title("Left trials")
% xticklabels([])
% xlabel("Linearized position (cm)")
ylabel(["Forward velocity"; "(cm/s)"])
% ylim([0,40])
% axis('square')

ax12 = subplot(2,2,2);
title("Right trials")
% xticklabels([])
% xlabel("Linearized Position (cm)")
% ylim([0,40])
yticklabels([])
% axis('square')

ax13 = subplot(2,2,3);
%title("Left Trials")
ylabel(["Angular velocity"; "(rad/s)"])
xlabel("Linearized position (cm)")
% ylim([-2,2])
% axis('square')
yline(0,'--','LineWidth',2);

ax14 = subplot(2,2,4);
%title("Right Trials")
xlabel("Linearized position (cm)")
% ylim([-2,2])
yticks([-1,0,1,2])
yticklabels([])
% axis('square')
yline(0,'--','LineWidth',2);

linkaxes([ax11,ax12])
linkaxes([ax13,ax14])

