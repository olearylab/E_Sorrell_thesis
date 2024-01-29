function [] = plot_error_cartoon(virmen_data,tbt_details,nbins,offsets)
% 05/05/2023

types_vec = [1,4,7,10];
colour_vec = [[0.4940 0.1840 0.5560];[0.4660 0.6740 0.1880]];
to_cm = 0.74;

% [error_mat,x_binned] = check_error_correction_07032023(virmen_data,tbt_details,nbins,offsets,false);
linearise_x = true;
x_vec = [5,6,7,15];
virmen_data(7,:) = wrapToPi(virmen_data(7,:));
[x_binned,centres] = bin_kin_data_05052023(virmen_data,x_vec,linearise_x,nbins);

mean_binned = nan.*ones(size(x_binned,2),size(x_binned,3),4);
for i = 1:4
    mean_binned(:,:,i) = squeeze(mean(x_binned(tbt_details(3,:)==types_vec(i),:,:),'omitnan'));
end

ITI = virmen_data(8,:);
cleaned_valid = clean_valid_data(ITI);
virmen_data = virmen_data(:,cleaned_valid);
trial_num = virmen_data(12,:);
xpos = to_cm*virmen_data(5,:);
ypos = to_cm*virmen_data(6,:);

figure

mazex = to_cm*[-0.1,0.1,0.1,1,1,15,15,-15,-15,-1,-1,-0.1,-0.1];
mazey = to_cm*[0,0,300,300,305,305,307.5,307.5,305,305,300,300,0];
plot(mazex,mazey,'Color',[0.5,0.5,0.5],'LineWidth',2)
hold on
axis equal
box off
axis off

% Example left trial
cur_trials = find(tbt_details(3,:) == types_vec(3));
n = 1;
plot(to_cm.*squeeze(x_binned(cur_trials(n),:,1)),to_cm.*squeeze(x_binned(cur_trials(n),:,2)),'Color',colour_vec(2,:),'LineWidth',2)

for b = 1:nbins
    quiver(to_cm.*squeeze(x_binned(cur_trials(n),b,1)),to_cm.*squeeze(x_binned(cur_trials(n),b,2)),(-sin(squeeze(mean_binned(b,3,1)))),(cos(squeeze(mean_binned(b,3,1)))),'Color',colour_vec(1,:),'LineWidth',3);
end

for b = 1:nbins
    quiver(to_cm.*squeeze(x_binned(cur_trials(n),b,1)),to_cm.*squeeze(x_binned(cur_trials(n),b,2)),(-sin(squeeze(x_binned(cur_trials(n),b,3)))),(cos(squeeze(x_binned(cur_trials(n),b,3)))),'Color',colour_vec(2,:),'LineWidth',3);
end