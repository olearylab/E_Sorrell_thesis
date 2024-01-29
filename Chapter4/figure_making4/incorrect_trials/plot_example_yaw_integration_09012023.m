function [] = plot_example_yaw_integration_09012023(virmen_data,t_num,start_dist,yaw_offset)
% 09/12/21

% Function for plotting single example incorrect BMI trial. Plot true ball
% yaw, true decoded view angle. End integrated yaw compared to view angle.

dt = 1/30;
ITI = virmen_data(8,:);
cleaned_valid = clean_valid_data(ITI);
% zdim = size(zdata,2);
trial_num_full = virmen_data(12,:);
% make sure ITI assigned to previous trial
trial_num_full = reassign_ITI_trial_num(cleaned_valid,trial_num_full);
trial_num_clean = trial_num_full(cleaned_valid);
num_trials = max(trial_num_clean);
virmen_data_clean = virmen_data(:,cleaned_valid);

circum = 64;
V = 0.32;
beta = 0.05*circum/V/2.5;
%offset = 1.495;
virmen_data_clean(15,:) = (virmen_data_clean(15,:)-yaw_offset)*-beta;

figure
subplot(2,2,2)
hold on
cur_trial = virmen_data_clean(:,trial_num_clean==t_num);
cur_time = (1:size(cur_trial,2))*dt;
start_index = find(cur_trial(6,:)>=start_dist,1);
plot(cur_time(1:start_index-1),cur_trial(7,1:start_index-1),'LineWidth',2,'Color',[0.5,0.5,0.5])
plot(cur_time(start_index:end),cur_trial(7,start_index:end),'LineWidth',2,'Color','k')
xline(cur_time(start_index),'LineWidth',2);
ylabel(["Decoded"; "View Angle (rad)"])
xticks([])

subplot(2,2,3)
hold on
plot(cur_time(1:start_index-1),cur_trial(15,1:start_index-1),'LineWidth',2,'Color',[0.5,0.5,0.5])
plot(cur_time(start_index:end),cur_trial(15,start_index:end),'LineWidth',2,'Color','b')
xline(cur_time(start_index),'LineWidth',2);
ylabel(["Ball Angular"; "Velocity (rad/s)"])
xlabel("Time (s)")

subplot(2,2,4)
hold on
plot(cur_time(1:start_index-1),cur_trial(7,1:start_index-1),'LineWidth',2,'Color',[0.5,0.5,0.5])
plot(cur_time(start_index:end),cur_trial(7,start_index:end),'LineWidth',2,'Color','k')
cur_va = cur_trial(7,start_index:end);
cur_yaw = (cur_trial(15,start_index:end)/(-beta)) + yaw_offset;
test_valid = ones(length(cur_yaw),1);
[va] = calc_va_alt_05_2021(cur_yaw,cur_va,yaw_offset,test_valid);
plot(cur_time(start_index:end),va,'LineWidth',2,'Color','b')
xline(cur_time(start_index),'LineWidth',2);
ylabel("View Angle (rad)")
xlabel("Time (s)")