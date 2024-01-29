function [stats_xc] = plot_neuron_correction_correlations(z_binned_cell,error_cell,tbt_cell,x_binned_cell,error_thresh)
% 27/06/2023

% Function for plotting correlation of neural activity with error
% correction

% Run on magnitude for now. Maybe consider looking at each direction
% separately or something.

% Maybe rerun using only tuned neurons or something?

% Just do not binned for now. Maybe move to binning by each relevant
% variable.

% Z should be normalised.
% Should bin by error in each position bin separately?
% Could do each neuron separately?

num_mice = size(error_cell,1);
num_days = size(error_cell,2);
types_vec = [1,4,7,10];
colour_vec = [[0.4940 0.1840 0.5560];[0.4660 0.6740 0.1880]];
n_error_bins = 20;
nbins = size(error_cell{1,1},2);

z_err_binned = nan.*ones(num_mice,num_days,n_error_bins,nbins,2);
all_err_centres = nan.*ones(n_error_bins,2,num_mice,num_days);

z_means1 = [];
errors_all1 = [];
signs_opposite1 = [];
x_all1 = [];
z_means2 = [];
errors_all2 = [];
signs_opposite2 = [];
x_all2 = [];

z_err_binned_av = nan.*ones(num_mice,num_days,n_error_bins,2);

for m = 1:num_mice
    for d = 1:num_days
        if ~isempty(error_cell{m,d})
            
            z_binned = z_binned_cell{m,d};
            error_mat = error_cell{m,d};
            tbt_details = tbt_cell{m,d};
            x_binned = x_binned_cell{m,d};
            
            % Just looking at bmi trials (both directions).
            for t = 1:2
                cur_trials = ismember(tbt_details(3,:),types_vec([1,2]+(t-1)*2));
                cur_error = error_mat(cur_trials,:);
                kept_bins = abs(cur_error)>error_thresh;
                cur_z = z_binned(cur_trials,:,:);
                cur_x = squeeze(x_binned(cur_trials,:,3)); % get just yaw
                signs_opposite = sign(cur_x) ~= sign(cur_error);
                % Could maybe make this consistent across all (just have to get
                % the max and min across all of them.
%                 [error_mag_bin,edges] = discretize(abs(cur_error),n_error_bins);
%                 all_err_centres(:,t,m,d) = (edges(2:end)+edges(1:end-1))/2;

%                 for i = 1:n_error_bins
%                     for j = 1:nbins
%                         z_err_binned(m,d,i,j,t) = mean(cur_z(error_mag_bin(:,j)==i,j,:),'all','omitnan');
%                     end
%                 end
                
%                 for i = 1:n_error_bins
%                     zz = permute(cur_z,[3,1,2]);
%                     zz = zz(:,:);
%                     z_err_binned_av(m,d,i,t) = mean(zz(:,error_mag_bin(:)==i),'all','omitnan');
%                 end
                
                pop_av = squeeze(mean(cur_z,3,'omitnan'));
                if t == 1
                    z_means1 = [z_means1;pop_av(kept_bins)];
                    errors_all1 = [errors_all1;abs(cur_error(kept_bins))];
                    x_all1 = [x_all1;abs(cur_x(kept_bins))];
                    signs_opposite1 = [signs_opposite1;signs_opposite(kept_bins)];
                else 
                    z_means2 = [z_means2;pop_av(kept_bins)];
                    errors_all2 = [errors_all2;abs(cur_error(kept_bins))];
                    x_all2 = [x_all2;abs(cur_x(kept_bins))];
                    signs_opposite2 = [signs_opposite2;signs_opposite(kept_bins)];
                end
            
            end
        end
    end
end
signs_opposite1 = signs_opposite1 == 1;
signs_opposite2 = signs_opposite2 == 1;

% z_err_binned_av = squeeze(mean(z_err_binned,4,'omitnan'));

%% Without error binning
% 
% errors_all1 = errors_all1.*(180/pi);
% errors_all2 = errors_all2.*(180/pi);
% 
% figure
% 
% % scatter(errors_all1,z_means1,'MarkerEdgeColor',colour_vec(1,:))
% scatter(errors_all2,z_means2,'MarkerEdgeColor',colour_vec(2,:))
% hold on
% 
% xlabel("View Angle Error (deg)")
% ylabel(["Mean BMI Trial";"Normalised Neural Activity"])
% yline(0,'--','LineWidth',2);


% [error_mag_bin,edges] = discretize(errors_all1,n_error_bins);
% centres = (edges(2:end)+edges(1:end-1))/2;
% final_binned = nan.*ones(n_error_bins,1);
% for i = 1:n_error_bins
%     final_binned(i) = mean(z_means1(error_mag_bin==i),'omitnan');
% end

% How to make it visible is tricky. Not sure how to do it with the same
% color.
% plot(centres(~isnan(final_binned)),final_binned(~isnan(final_binned)),'-o','Color','k','LineWidth',3)
% 
% [error_mag_bin,edges] = discretize(errors_all2,n_error_bins);
% centres = (edges(2:end)+edges(1:end-1))/2;
% final_binned = nan.*ones(n_error_bins,1);
% for i = 1:n_error_bins
%     final_binned(i) = mean(z_means2(error_mag_bin==i),'omitnan');
% end
% 
% plot(centres(~isnan(final_binned)),final_binned(~isnan(final_binned)),'--o','Color','k','LineWidth',3)

%% Other stats

% % Spearman
% [stats.s1,stats.sp1] = corr(errors_all1,z_means1,'Type','Spearman','rows','pairwise');
% [stats.s2,stats.sp2] = corr(errors_all2,z_means2,'Type','Spearman','rows','pairwise');
% 
% % [stats.sb1,stats.sbp1] = corr(er_all(g_all==1),z_all(g_all==1),'Type','Spearman','rows','pairwise');
% % [stats.sb2,stats.sbp2] = corr(er_all(g_all==2),z_all(g_all==2),'Type','Spearman','rows','pairwise');
% 
% % Pearson
% [stats.rho1,stats.p1] = corr(errors_all1,z_means1,'rows','pairwise');
% [stats.rho2,stats.p2] = corr(errors_all2,z_means2,'rows','pairwise');
% 
% % [stats.rhob1,stats.pb1] = corr(er_all(g_all==1),z_all(g_all==1),'rows','pairwise');
% % [stats.rhob2,stats.pb2] = corr(er_all(g_all==2),z_all(g_all==2),'rows','pairwise');
% 
% % Linear regression - careful

%% z with absolute movements
% x_all1 = x_all1.*(180/pi);
x_all2 = x_all2.*(180/pi);

% figure
% 
% scatter(x_all1,z_means1,'MarkerEdgeColor',colour_vec(1,:))
% hold on
% scatter(x_all2,z_means2,'MarkerEdgeColor',colour_vec(2,:))
% 
% xlabel("View Angle Velocity (deg/s)")
% ylabel(["Mean BMI Trial";"Normalised Neural Activity"])
% yline(0,'--','LineWidth',2);
% 
% 
% [x_mag_bin,edges] = discretize(x_all1,n_error_bins);
% centres = (edges(2:end)+edges(1:end-1))/2;
% final_binned = nan.*ones(n_error_bins,1);
% for i = 1:n_error_bins
%     final_binned(i) = mean(z_means1(x_mag_bin==i),'omitnan');
% end
% 
% % How to make it visible is tricky. Not sure how to do it with the same
% % color.
% plot(centres(~isnan(final_binned)),final_binned(~isnan(final_binned)),'-o','Color','k','LineWidth',3)
% 
% [x_mag_bin,edges] = discretize(x_all2,n_error_bins);
% centres = (edges(2:end)+edges(1:end-1))/2;
% final_binned = nan.*ones(n_error_bins,1);
% for i = 1:n_error_bins
%     final_binned(i) = mean(z_means2(x_mag_bin==i),'omitnan');
% end
% 
% plot(centres(~isnan(final_binned)),final_binned(~isnan(final_binned)),'--o','Color','k','LineWidth',3)
% 
% %% Run ANCOVA
% 
% x_all = [x_all1;x_all2];
% z_means = [z_means1;z_means2];
% groups_vec = [ones(length(x_all1),1);2.*ones(length(x_all2),1)];
% 
% aoctool(x_all,z_means,groups_vec);
% %% Other stats
% 
% % Spearman
% [stats_x.s1,stats_x.sp1] = corr(x_all1,z_means1,'Type','Spearman','rows','pairwise');
% [stats_x.s2,stats_x.sp2] = corr(x_all2,z_means2,'Type','Spearman','rows','pairwise');
% 
% % Pearson
% [stats_x.rho1,stats_x.p1] = corr(x_all1,z_means1,'rows','pairwise');
% [stats_x.rho2,stats_x.p2] = corr(x_all2,z_means2,'rows','pairwise');
% 

%% z with error correcting movements
figure

% scatter(x_all1(signs_matching1),z_means1(signs_matching1),'MarkerEdgeColor',colour_vec(1,:))
scatter(x_all2(signs_opposite2),z_means2(signs_opposite2),'.','MarkerEdgeColor',colour_vec(2,:))
hold on

xlabel("View Angle Velocity Magnitude (deg/s)")
ylabel(["Mean BMI Trial";"Normalised Neural Activity"])
yline(0,'--','LineWidth',2);

mdl = fitlm(x_all2(signs_opposite2),z_means2(signs_opposite2));
xx = linspace(min(x_all2(signs_opposite2)),max(x_all2(signs_opposite2)),20)';
[ypred,yci] = predict(mdl,xx);
plot(xx,yci(:,1),'--','LineWidth',2,'Color','k')
plot(xx,yci(:,2),'--','LineWidth',2,'Color','k')
plot(xx,ypred,'-','LineWidth',2,'Color','k')
h = fill([xx',fliplr(xx')],[yci(:,2)',fliplr(yci(:,1)')],'k','LineStyle','none');
set(h,'facealpha',.1)

box off
axis('square')
title(["Population neural activity";"vs heading correcting movement"]);


% [xc_mag_bin,edges] = discretize(x_all1(signs_matching1),n_error_bins);
% centres = (edges(2:end)+edges(1:end-1))/2;
% final_binned = nan.*ones(n_error_bins,1);
% zz = z_means1(signs_matching1);
% for i = 1:n_error_bins
%     final_binned(i) = mean(zz(xc_mag_bin==i),'omitnan');
% end
% 
% % How to make it visible is tricky. Not sure how to do it with the same
% % color.
% plot(centres(~isnan(final_binned)),final_binned(~isnan(final_binned)),'-o','Color','k','LineWidth',3)
% 
% [xc_mag_bin,edges] = discretize(x_all2(signs_matching2),n_error_bins);
% centres = (edges(2:end)+edges(1:end-1))/2;
% final_binned = nan.*ones(n_error_bins,1);
% zz = z_means2(signs_matching2);
% for i = 1:n_error_bins
%     final_binned(i) = mean(zz(xc_mag_bin==i),'omitnan');
% end
% 
% plot(centres(~isnan(final_binned)),final_binned(~isnan(final_binned)),'--o','Color','k','LineWidth',3)

%% Run ANCOVA
% 
% xc_all = [x_all1(signs_matching1);x_all2(signs_matching2)];
% zc_means = [z_means1(signs_matching1);z_means2(signs_matching2)];
% groups_vec = [ones(length(x_all1(signs_matching1)),1);2.*ones(length(x_all2(signs_matching2)),1)];
% 
% aoctool(xc_all,zc_means,groups_vec);
%% Other stats

% Spearman
[stats_xc.s1,stats_xc.sp1] = corr(x_all1(signs_opposite1),z_means1(signs_opposite1),'Type','Spearman','rows','pairwise');
[stats_xc.s2,stats_xc.sp2] = corr(x_all2(signs_opposite2),z_means2(signs_opposite2),'Type','Spearman','rows','pairwise');

% Pearson
[stats_xc.rho1,stats_xc.p1] = corr(x_all1(signs_opposite1),z_means1(signs_opposite1),'rows','pairwise');
[stats_xc.rho2,stats_xc.p2] = corr(x_all2(signs_opposite2),z_means2(signs_opposite2),'rows','pairwise');
