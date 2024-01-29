function [error_cell,x_cell,ro_all,sign_opposite,h_boots] = run_many_check_error_correction_norm(virmen_cell,tbt_cell,virmen_train_cell,tbt_train_cell,nbins,error_thresh)
% 09/09/2023

num_mice = size(virmen_cell,1);
num_days = size(virmen_cell,2);

types_vec = [1,4,7,10];
% types_vec = [2,5,8,11];
num_shuffles = 100;

offsets = [1.4953,1.4961;1.4953,1.4961;1.4953,1.4961;1.4804,1.4844;1.4804,1.4844];

error_cell = cell(num_mice,num_days);
x_cell = cell(num_mice,num_days);

for m = 1:num_mice
    [mean_binned,std_binned] = calculate_mean_ball_va(virmen_train_cell{m},tbt_train_cell{m},nbins,offsets(m,:));
    for d = 1:num_days
        if ~isempty(virmen_cell{m,d})
 
            % [error_cell{m,d},x_cell{m,d}] = check_error_correction_07032023(virmen_cell{m,d},tbt_cell{m,d},nbins,offsets(m,:),false);
            [error_cell{m,d},x_cell{m,d}] = check_error_correction_norm(virmen_cell{m,d},tbt_cell{m,d},mean_binned,std_binned,nbins,offsets(m,:),false);
            
        end
    end
end

%% Plotting
% Combine left and right
plot_start = 2; % 0 for ball trials, 2 for BMI trials
n_error_bins = 6;
ro_all = nan.*ones(num_mice,num_days);
sign_opposite = nan.*ones(num_mice,num_days);
sign_opposite_bins = nan.*ones(nbins,num_mice,num_days);
ro_bins = nan.*ones(nbins,num_mice,num_days);
ro_bins_shuff = nan.*ones(num_shuffles,nbins,num_mice,num_days);
sign_opposite_abs_error_bins = nan.*ones(n_error_bins,num_mice,num_days);
num_samp_abs_error_bins = nan.*ones(n_error_bins,num_mice,num_days);
abs_error_centres = nan.*ones(n_error_bins,num_mice,num_days);

err_range = [1,3];
sign_opposite_bins_range = nan.*ones(nbins,num_mice,num_days);

sign_opposite_error_bins = nan.*ones(n_error_bins,num_mice,num_days);
num_samp_error_bins = nan.*ones(n_error_bins,num_mice,num_days);
error_centres = nan.*ones(n_error_bins,num_mice,num_days);
for m = 1:num_mice
    for d = 1:num_days
        if ~isempty(virmen_cell{m,d})
            
            error_mat = error_cell{m,d};
            x_binned = x_cell{m,d};
            tbt_details = tbt_cell{m,d};
            
            % figure

            cur_trials = find(ismember(tbt_details(3,:),types_vec([plot_start+1,plot_start+2])));
            % cur_trials = find(ismember(tbt_details(3,:),[7,8,10,11])); % For both correct and incorrect
            cur_errors = error_mat(cur_trials,:);
            cur_yaw = squeeze(x_binned(cur_trials,:,3));
            ro = corr(cur_errors(abs(cur_errors)>error_thresh),cur_yaw(abs(cur_errors)>error_thresh));
            ro_all(m,d) = ro;
            for b = 1:nbins
                ro_bins(b,m,d) = corr(cur_errors(abs(cur_errors(:,b))>error_thresh,b),cur_yaw(abs(cur_errors(:,b))>error_thresh,b));
            end
            
            % Calculate shuffled correlations
            for n = 1:num_shuffles
                cur_yaw_shuff = cur_yaw(:);
                cur_yaw_shuff = cur_yaw_shuff(randperm(length(cur_yaw_shuff)));
                cur_yaw_shuff = reshape(cur_yaw_shuff,size(cur_yaw,1),size(cur_yaw,2));
                for b = 1:nbins
                    ro_bins_shuff(n,b,m,d) = corr(cur_errors(abs(cur_errors(:,b))>error_thresh,b),cur_yaw_shuff(abs(cur_errors(:,b))>error_thresh,b));
                end
            end
                
            
            % Convert to degrees - not any more since normalised
            % cur_errors = cur_errors.*(180/pi);
            
            % Check matching sign
            error_sign = cur_errors>0;
            yaw_sign = cur_yaw>0;
            sign_opposite(m,d) = sum(error_sign(abs(cur_errors)>error_thresh)~=yaw_sign(abs(cur_errors)>error_thresh))/length(error_sign(abs(cur_errors)>error_thresh));
            % Check signs bin-wise
            % Check matching sign in specific error range
            in_range = (abs(cur_errors)>err_range(1))&(abs(cur_errors)<err_range(2));
            for b = 1:nbins
                sign_opposite_bins(b,m,d) = sum(error_sign(abs(cur_errors(:,b))>error_thresh,b)~=yaw_sign(abs(cur_errors(:,b))>error_thresh,b))./size(error_sign(abs(cur_errors(:,b))>error_thresh,b),1);
                sign_opposite_bins_range(b,m,d) = sum(error_sign(in_range(:,b),b)~=yaw_sign(in_range(:,b),b))./size(error_sign(in_range(:,b),b),1); 
            end
            
            % Check matching sign binned by magnitude of error

            [error_mag_bin,edges] = discretize(abs(cur_errors),n_error_bins);
            abs_error_centres(:,m,d) = (edges(2:end)+edges(1:end-1))/2;
            sign_opposite_all = error_sign~=yaw_sign;
            for b = 1:n_error_bins
                sign_opposite_abs_error_bins(b,m,d) = sum(sign_opposite_all(error_mag_bin==b))./sum(error_mag_bin(:)==b);
                num_samp_abs_error_bins(b,m,d) = sum(error_mag_bin(:)==b);
            end
                
            % Check matching sign binned by error

            [error_mag_bin,edges] = discretize(cur_errors,n_error_bins);
            error_centres(:,m,d) = (edges(2:end)+edges(1:end-1))/2;
            sign_opposite_all = error_sign~=yaw_sign;
            for b = 1:n_error_bins
                sign_opposite_error_bins(b,m,d) = sum(sign_opposite_all(error_mag_bin==b))./sum(error_mag_bin(:)==b);
                num_samp_error_bins(b,m,d) = sum(error_mag_bin(:)==b);
            end
        end
            
    end
end

%% Get centres
[x_binned,centres] = bin_kin_data(virmen_cell{1,1},6,true,nbins);
centres = centres*0.74;

%% PLot
colour_vec = [[0.4940 0.1840 0.5560];[0.4660 0.6740 0.1880]];
% figure
% violin(ro_all(:),'facecolor',colour_vec(2,:));
% hold on
% % scatter(ones(num_mice*num_days,1),ro_all(:),100,'filled','k')
% % xlim([0.5,2.5])
% ylim([0,1])
% ylabel("Correlation")
% % xlabel("Trial Type")
% xticks([1])
% % xticklabels(["BMI Left";"BMI Right"])
% xticklabels(["BMI Trials"])
% title(["Comparison of ball angular velocity"; "and heading correction directions"])
% % title(["Correlation between"; "ball angular velocity and heading correction vector"])
% axis('square')
% box off

% figure
% violin(sign_matching(:),'facecolor',colour_vec(2,:));
% hold on
% % scatter(ones(num_mice*num_days,1),sign_matching(1,:),100,'filled','k')
% 
% % yline(0.5,'--','LineWidth',2);
% % xlim([0.5,2.5])
% ylim([0.5,1])
% title(["Comparison of ball angular velocity"; "and heading correction directions"])
% ylabel(["Fraction matching signs"])
% % xlabel("Trial Type")
% xticks([1])
% xticklabels(["BMI Trials"])
% % title(["Comparison of"; "yaw and error direction"])
% axis('square')
% box off

% figure
% % titles = ["BMI Left";"BMI Right"];
% for n = 1:nbins
%     plot(centres(n).*ones(num_mice*num_days,1),squeeze(sign_matching_bins(n,:)),'o','Color','k')
%     hold on
% end
% plot(centres,squeeze(mean(sign_matching_bins(:,:),2,'omitnan')),'LineWidth',2,'Color','k')
% %xlim([0,nbins+1])
% ylim([0,1])
% yline(0.5,'--','LineWidth',2);
% title(["Comparison of ball angular velocity"; "and heading correction directions"])
% xlabel("Linearised Position (cm)")
% ylabel(["Fraction matching signs"])
% axis('square')
% box off

% Alternative with lines
% figure
% % titles = ["BMI Left";"BMI Right"];
% sign_matching_bins = sign_matching_bins(:,:);
% for n = 1:size(sign_matching_bins,2)
%     plot(centres,squeeze(sign_matching_bins(:,n)),'Color',[0.5,0.5,0.5])
%     hold on
% end
% plot(centres,squeeze(mean(sign_matching_bins(:,:),2,'omitnan')),'LineWidth',2,'Color','k')
% %xlim([0,nbins+1])
% ylim([0,1])
% yline(0.5,'--','LineWidth',2);
% title(["Comparison of ball angular velocity"; "and heading deviation"])
% xlabel("Linearised Position (cm)")
% ylabel(["Fraction opposite signs"])
% axis('square')
% box off
% Final figure of above probably just wants means and standard deviations or
% convidence interval, or bootstrap confidence interval?

% Alternative with correlations
% figure
% ro_bins = ro_bins(:,:);
% for n = 1:size(ro_bins,2)
%     plot(centres,squeeze(ro_bins(:,n)),'Color',[0.5,0.5,0.5])
%     hold on
% end
% plot(centres,squeeze(mean(ro_bins(:,:),2,'omitnan')),'LineWidth',2,'Color','k')
% %xlim([0,nbins+1])
% ylim([-1,1])
% yline(0,'--','LineWidth',2);
% title(["Correlation between heading deviation"; "and ball angular velocity"])
% xlabel("Linearised Position (cm)")
% ylabel(["Pearson Correlation"])
% axis('square')
% box off

%% Alternative with correlations and bootstrapping
% Maybe add in shuffled results to run stats test on?

boot_samps = 1000;
num_trials = 4;
% Calculates probability that data2 > data1. 
% [p_boot, bootstats, bootstats_center, bootstats_sem] = get_bootstrap_results_equalsamples(squeeze(p_R2_all_mat(:,:,2)),squeeze(p_R2_all_mat(:,:,4)),boot_samps,num_trials,'mean');
orig_ro_bins = ro_bins;
for m = 1:num_mice
    d_ind = 0;
    for d = 1:num_days
        if ~isnan(orig_ro_bins(1,m,d)) 
            d_ind = d_ind+1;
            ro_bins(:,m,d_ind) = orig_ro_bins(:,m,d);

        end
    end
    if d_ind<num_days
        ro_bins(:,m,d_ind+1:num_days) = nan;
    end
end

orig_ro_bins_shuff = ro_bins_shuff;
for m = 1:num_mice
    d_ind = 0;
    for d = 1:num_days
        if ~isnan(orig_ro_bins_shuff(1,1,m,d)) 
            d_ind = d_ind+1;
            ro_bins_shuff(:,:,m,d_ind) = orig_ro_bins_shuff(:,:,m,d);

        end
    end
    if d_ind<num_days
        ro_bins_shuff(:,:,m,d_ind+1:num_days) = nan;
    end
end

ro_bins_shuff = squeeze(mean(ro_bins_shuff));

all_centres = NaN(nbins,2);
all_sems = NaN(nbins,2);
all_p_boot = NaN(nbins,1);
for b = 1:nbins
    [all_p_boot(b),all_centres(b,:),all_sems(b,:)] = run_H_boot_ets(squeeze(ro_bins(b,:,:)), squeeze(ro_bins_shuff(b,:,:)),true);
end

h_boots.all_p_boot = all_p_boot;
h_boots.all_centres = all_centres;
h_boots.all_sems = all_sems;
% 
% bootstats_s = NaN(nbins,boot_samps);
% for b = 1:nbins
%     bootstats_s(b,:) = get_bootstrapped_equalsamples(squeeze(mean(ro_bins_shuff(:,b,:,:))),boot_samps,num_trials,'mean');
% end
% 
% bootstats = NaN(nbins,boot_samps);
% for b = 1:nbins
%     bootstats(b,:) = get_bootstrapped_equalsamples(squeeze(ro_bins(b,:,:)),boot_samps,num_trials,'mean');
% end

% CI_vals = [0.5,99.5]; % [2.5,97.5];
% CI_vals = [2.5,97.5];
% 
% % CIs, upper/lower x ntypes x nbins x zdim
% CIs_all = zeros(2,nbins);
% for i = 1:nbins
%     CIs_all(:,i) = prctile(bootstats(i,:),CI_vals);
% end

figure

% boostrap stds instead
boot_stds = all_sems;
lims_all = zeros(2,nbins,2);
lims_all(1,:,:) = all_centres - boot_stds;
lims_all(2,:,:) = all_centres + boot_stds;

% h = fill([centres,fliplr(centres)],[CIs_all(1,:),fliplr(CIs_all(2,:))],'k','EdgeColor','none');
h = fill([centres,fliplr(centres)],[squeeze(lims_all(1,:,1)),fliplr(squeeze(lims_all(2,:,1)))],'k','EdgeColor','none');
set(h,'facealpha',.3)
hold on

plot(centres,all_centres(:,1),'LineWidth',2,'Color','k')


% CI_vals = [0.5,99.5]; % [2.5,97.5];
% CI_vals = [2.5,97.5];
% 
% % CIs, upper/lower x ntypes x nbins x zdim
% CIs_all_s = zeros(2,nbins);
% for i = 1:nbins
%     CIs_all_s(:,i) = prctile(bootstats_s(i,:),CI_vals);
% end

% boostrap stds instead
% boot_stds_s = std(bootstats_s,0,2);
% lims_all_s = zeros(2,nbins);
% lims_all_s(1,:) = mean(bootstats_s,2,'omitnan') - boot_stds_s;
% lims_all_s(2,:) = mean(bootstats_s,2,'omitnan') + boot_stds_s;

% h = fill([centres,fliplr(centres)],[CIs_all_s(1,:),fliplr(CIs_all_s(2,:))],'k','EdgeColor','none');
h = fill([centres,fliplr(centres)],[squeeze(lims_all(1,:,2)),fliplr(squeeze(lims_all(2,:,2)))],'k','EdgeColor','none');
set(h,'facealpha',.3)
hold on

plot(centres,all_centres(:,2),'--','LineWidth',2,'Color','k')

%xlim([0,nbins+1])
ylim([-1,1])
% yline(0,'--','LineWidth',2);
title(["Correlation between heading deviation"; "and ball angular velocity"])
xlabel("Linearised Position (cm)")
ylabel(["Pearson Correlation"])
axis('square')
box off
%%

% In range
% figure
% % titles = ["BMI Left";"BMI Right"];
% sign_matching_bins_range = sign_matching_bins_range(:,:);
% for n = 1:size(sign_matching_bins_range,2)
%     plot(centres,squeeze(sign_matching_bins_range(:,n)),'Color',[0.5,0.5,0.5])
%     hold on
% end
% plot(centres,squeeze(mean(sign_matching_bins_range(:,:),2,'omitnan')),'LineWidth',2,'Color','k')
% %xlim([0,nbins+1])
% ylim([0,1])
% yline(0.5,'--','LineWidth',2);
% title(["Comparison of ball angular velocity"; "and heading deviation";"25<Error Magnitude<50"])
% xlabel("Linearised Position (cm)")
% ylabel(["Fraction opposite signs"])
% axis('square')
% box off
% 
% % abs error bin plot
% sign_matching_abs_error_bins = sign_matching_abs_error_bins(:,:);
% abs_error_centres = abs_error_centres(:,:);
% figure
% % titles = ["BMI Left";"BMI Right"];
% for n = 1:size(sign_matching_abs_error_bins,2)
%     plot(squeeze(abs_error_centres(:,n)),squeeze(sign_matching_abs_error_bins(:,n)),'-o','Color',[0.5,0.5,0.5])
%     hold on
% end
% 
% % Get binned means
% [e_disc,edges] = discretize(abs_error_centres,n_error_bins);
% mean_binned = nan.*ones(n_error_bins,1);
% eb_centres = (edges(2:end)+edges(1:end-1))/2;
% for b = 1:n_error_bins
%     mean_binned(b) = mean(sign_matching_abs_error_bins(e_disc==b),'omitnan');
% end
% 
% plot(eb_centres,mean_binned,'-o','LineWidth',2,'Color','k')
% %xlim([0,nbins+1])
% % ylim([0,1])
% yline(0.5,'--','LineWidth',2);
% title(["Comparison of ball angular velocity"; "and heading correction directions"])
% xlabel("Magnitude of Heading Deviation (deg)")
% ylabel(["Fraction matching signs"])
% axis('square')
% box off
% 
% % error bin plot
% sign_matching_error_bins = sign_matching_error_bins(:,:);
% error_centres = error_centres(:,:);
% figure
% % titles = ["BMI Left";"BMI Right"];
% 
% for n = 1:size(sign_matching_error_bins,2)
%     plot(squeeze(error_centres(:,n)),squeeze(sign_matching_error_bins(:,n)),'-o','Color',[0.5,0.5,0.5])
%     hold on
% end
% 
% % Get binned means
% [e_disc,edges] = discretize(error_centres,n_error_bins+1);
% mean_binned = nan.*ones(n_error_bins+1,1);
% eb_centres = (edges(2:end)+edges(1:end-1))/2;
% for b = 1:n_error_bins+1
%     mean_binned(b) = mean(sign_matching_error_bins(e_disc==b),'omitnan');
% end
% 
% plot(eb_centres,mean_binned,'-o','LineWidth',2,'Color','k')
% %xlim([0,nbins+1])
% % ylim([0,1])
% yline(0.5,'--','LineWidth',2);
% xline(0,'--','LineWidth',2);
% title(["Comparison of ball angular velocity"; "and heading correction directions"])
% xlabel("Heading Deviation (deg)")
% ylabel(["Fraction matching signs"])
% axis('square')
% box off

%% get p_boots
% all_p_boot = nan.*ones(nbins,1);
% for b = 1:nbins
%     all_p_boot(b) = get_direct_prob(bootstats_s(b,:),bootstats(b,:));
% end