function [] = plot_example_neuron_correction_correlation(z_binned_cell,x_binned_cell,error_cell,tbt_cell,z_error_corrs_cell,ex_md)
% 27/06/2023

% plot scatter of correlation for a single neuron between activity and
% error correcting behaviour.
% Could plot as "tuning curve" line instead of scatter?
% Would have to bin by error to plot as error tuning curve
types_vec = [1,4,7,10];

z_binned = z_binned_cell{ex_md(1),ex_md(2)};
error_mat = error_cell{ex_md(1),ex_md(2)};
tbt_details = tbt_cell{ex_md(1),ex_md(2)};
z_error_corrs = z_error_corrs_cell{ex_md(1),ex_md(2)};
x_binned = x_binned_cell{ex_md(1),ex_md(2)};

[~,n_ind] = max(z_error_corrs(:,3));

cur_trials = ismember(tbt_details(3,:),types_vec([3,4]));

zz = squeeze(z_binned(cur_trials,:,n_ind));
xx = squeeze(x_binned(cur_trials,:,3)).*(180/pi);
ee = error_mat(cur_trials,:).*180/pi;
signs_matching = sign(ee) == sign(xx);

colour_vec = [[0.4940 0.1840 0.5560];[0.4660 0.6740 0.1880]];

figure
scatter(abs(xx(signs_matching)),zz(signs_matching),'filled','MarkerEdgeColor',colour_vec(2,:),'MarkerFaceColor',colour_vec(2,:))

hold on

mdl = fitlm(abs(xx(signs_matching)),zz(signs_matching));
xx = linspace(min(abs(xx(signs_matching))),max(abs(xx(signs_matching))),20)';
[ypred,yci] = predict(mdl,xx);
plot(xx,yci(:,1),'--','LineWidth',2,'Color','k')
plot(xx,yci(:,2),'--','LineWidth',2,'Color','k')
plot(xx,ypred,'-','LineWidth',2,'Color','k')
h = fill([xx',fliplr(xx')],[yci(:,2)',fliplr(yci(:,1)')],'k','LineStyle','none');
set(h,'facealpha',.1)

box off
xlabel("View Angle Velocity Magnitude (deg/s)")
ylabel(["Mean BMI Trial";"Normalised Neural Activity"])
yline(0,'--','LineWidth',2);
title(["Example Neural Activity";"vs Heading Correcting Movements"])
axis('square')