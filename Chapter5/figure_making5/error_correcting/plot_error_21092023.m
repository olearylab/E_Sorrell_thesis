function [h_boots] = plot_error_21092023(error_cell,tbt_cell,nbins,centres)
% 21/09/2023

% Plot view angle "error" as a function of position

% Should me error_cell_norm

num_mice = size(error_cell,1);
num_days = size(error_cell,2);
types_vec = [1,4,7,10];

centres = centres*0.74;
offsets_mat = [1.4953,1.4961;1.4953,1.4961;1.4953,1.4961;1.4804,1.4844;1.4804,1.4844];
figure

mean_errors = nan.*ones(nbins,num_mice,num_days);

for m = 1:num_mice
    for d = 1:num_days
        if ~isempty(error_cell{m,d})
            
            tbt_details = tbt_cell{m,d};
            error_mat = error_cell{m,d};
            cur_trials = ismember(tbt_details(3,:),types_vec([3,4]));
            mean_errors(:,m,d) = mean(abs(error_mat(cur_trials,:)));
            
        end
    end
end

figure

mean_errors = mean_errors(:,:);
for n = 1:size(mean_errors,2)
    plot(centres,squeeze(mean_errors(:,n)),'Color',[0.5,0.5,0.5])
    hold on
end
plot(centres,squeeze(mean(mean_errors(:,:),2,'omitnan')),'LineWidth',2,'Color','k')
%xlim([0,nbins+1])
% ylim([0,1])
% yline(0.5,'--','LineWidth',2);
title("View Angle Errors")
xlabel("Linearised Position (cm)")
ylabel(["Error Magnitude (deg)"])
axis('square')
box off
 
%% H Bootstrap

all_centres = NaN(nbins,2);
all_sems = NaN(nbins,2);
all_p_boot = NaN(nbins,1);
for b = 1:nbins
    [all_p_boot(b),all_centres(b,:),all_sems(b,:)] = run_H_boot_ets(squeeze(mean_errors(b,:,:)), squeeze(mean_errors(b,:,:)),true);
end

h_boots.all_p_boot = all_p_boot;
h_boots.all_centres = all_centres;
h_boots.all_sems = all_sems;

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

%xlim([0,nbins+1])
% ylim([-1,1])
% yline(0,'--','LineWidth',2);
title(["Heading deviations"; "along the T-maze"])
xlabel("Linearised position (cm)")
ylabel(["Mean heading deviation" "magnitude (a.u.)"])
axis('square')
box off
% xticks([0,200])
ylim([0,3])
yticks([0,1,2,3])