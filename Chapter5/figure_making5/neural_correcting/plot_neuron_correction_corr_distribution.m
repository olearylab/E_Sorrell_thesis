function [n_ind] = plot_neuron_correction_corr_distribution(z_binned_cell,x_binned_cell,error_cell,tbt_cell,z_correction_corrs_cell,ex_md,error_thresh)
% 27/06/2023

% plot scatter of correlation for a single neuron between activity and
% error correcting behaviour.
% Could plot as "tuning curve" line instead of scatter?
% Would have to bin by error to plot as error tuning curve

% With distributions of neuron correlaitons.
% Need z_error_corrs_cell to be correlation with heading correcting
% movements, not heading deviations.

types_vec = [1,4,7,10];

num_mice = size(z_binned_cell,1);
num_days = size(z_binned_cell,2);

z_binned = z_binned_cell{ex_md(1),ex_md(2)};
error_mat = error_cell{ex_md(1),ex_md(2)};
tbt_details = tbt_cell{ex_md(1),ex_md(2)};
z_correction_corrs = z_correction_corrs_cell{ex_md(1),ex_md(2)};
x_binned = x_binned_cell{ex_md(1),ex_md(2)};

% only one column (didn't bother with directions separate)
% [~,n_ind] = max(z_correction_corrs(:,3));
[~,n_ind] = max(z_correction_corrs);

cur_trials = ismember(tbt_details(3,:),types_vec([3,4]));

zz = squeeze(z_binned(cur_trials,:,n_ind));
xx = squeeze(x_binned(cur_trials,:,3)).*(180/pi);
ee = error_mat(cur_trials,:); % .*180/pi;
signs_opposite = sign(ee) ~= sign(xx);

kept_bins = abs(ee)>error_thresh;

zz = zz(kept_bins);
xx = xx(kept_bins);
signs_opposite = signs_opposite(kept_bins);

colour_vec = [[0.4940 0.1840 0.5560];[0.4660 0.6740 0.1880]];

figure
scatter(abs(xx(signs_opposite)),zz(signs_opposite),'filled','MarkerEdgeColor',colour_vec(2,:),'MarkerFaceColor',colour_vec(2,:))

hold on

mdl = fitlm(abs(xx(signs_opposite)),zz(signs_opposite));
xx = linspace(min(abs(xx(signs_opposite))),max(abs(xx(signs_opposite))),20)';
[ypred,yci] = predict(mdl,xx);
plot(xx,yci(:,1),'--','LineWidth',2,'Color','k')
plot(xx,yci(:,2),'--','LineWidth',2,'Color','k')
plot(xx,ypred,'-','LineWidth',2,'Color','k')
h = fill([xx',fliplr(xx')],[yci(:,2)',fliplr(yci(:,1)')],'k','LineStyle','none');
set(h,'facealpha',.1)

box off
xlabel("View angle velocity magnitude (deg/s)")
ylabel(["Mean BMI trial";"normalised neural activity"])
yline(0,'--','LineWidth',2);
title(["Example neural activity";"vs heading correcting movements"])
axis('square')

%% Distributions 
num_sess = num_mice*num_days;
v_all = [];
v_sess = cell(1,num_sess-1);
sess_ind = 1;
for m = 1:num_mice
    for d = 1:num_days
        if ~isempty(z_correction_corrs_cell{m,d})
            
            z_correction_corrs = z_correction_corrs_cell{m,d};
            v_all = [v_all;z_correction_corrs];
            v_sess{sess_ind} = z_correction_corrs;
            sess_ind = sess_ind + 1;
        end
    end
end

figure
violin(v_all,'facecolor',colour_vec(2,:),'medc',[],'edgecolor',[]);
yline(0,'--','LineWidth',2);
title(["Distribution of Activity Correlations";"with Corrective Movements"])
ylabel("Pearson Correlations")
xticks([1])
box off
axis('square')


figure
violin(v_sess,'facecolor',colour_vec(2,:),'medc',[],'edgecolor',[]);
title(["Distribution of Activity Correlations";"with Corrective Movements"])
ylabel("Pearson Correlations")
xlabel("Session")

yline(0,'--','LineWidth',2);
box off
axis('square')