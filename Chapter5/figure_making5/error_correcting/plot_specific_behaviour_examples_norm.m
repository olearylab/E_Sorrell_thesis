function [] = plot_specific_behaviour_examples_norm(virmen_data,tbt_details,nbins,offsets,mean_binned,std_binned,centres)
% 01/05/2023

% Plot mean ball trial trajectory, and BMI trial trajectory for cartoon of
% how error is calcualted.

plot_res = false;
types_vec = [1,4,7,10];
centres = centres*0.74;

[error_mat,x_binned] = BMI_correct_errors_norm(virmen_data,tbt_details,mean_binned,std_binned,nbins,offsets);

ex_x_binned = x_binned(tbt_details(3,:)==types_vec(3),:,:);
ex_error_mat = error_mat(tbt_details(3,:)==types_vec(3),:);

ex_ind = 6;
figure
yyaxis left
colour_vec = [[0.4940 0.1840 0.5560];[0.4660 0.6740 0.1880]];
plot(centres,ex_error_mat(ex_ind,:),'Color','k','LineWidth',3)
ylabel("Normalised heading deviation (a.u.)")
ylim([-5,5])
hold on
yyaxis right
plot(centres,ex_x_binned(ex_ind,:,3),'-','Color',colour_vec(2,:),'LineWidth',3)
hold on
plot(centres,ex_x_binned(ex_ind,:,5),':','Color',colour_vec(2,:),'LineWidth',3)
ylabel("Angular velocity (rad/s)")

xlabel("Linearised position (cm)")
title("Heading correction example")
yline(0,'--','LineWidth',2);
ylim([-1.5,1.5])
% legend("Average of Ball Left Trials","Example BMI Left Trial")
box off
axis('square')
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = colour_vec(2,:);

% Plot many trials for current mouse and session
% (need different funciton to do many mice and sessions at once.
figure
num_ex = size(ex_x_binned,1);
num_plot = min([num_ex,9]);
% look at 9 examples for now (if there are that many)
for i = 1:num_plot
    subplot(3,3,i)
    yyaxis left
    plot(centres,ex_error_mat(i,:),'Color','k','LineWidth',3)
    % ylabel("View Angle Error (rad)")
    ylim([-5,5])
    hold on
    yyaxis right
    plot(centres,ex_x_binned(i,:,3),'Color',colour_vec(2,:),'LineWidth',3)

    % ylabel("Angular Velocity (rad/s)")

    % xlabel("Linearised Position (cm)")
    % title("View Angle Error Correction Example")
    yline(0,'--','LineWidth',2);
    ylim([-1.6,1.6])
    % legend("Average of Ball Left Trials","Example BMI Left Trial")
    box off
    ax = gca;
    ax.YAxis(1).Color = 'k';
    ax.YAxis(2).Color = colour_vec(2,:);
end
