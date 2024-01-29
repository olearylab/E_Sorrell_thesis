function [] = plot_early_fraction_correct(tbt_cell,av_num)
% 09/06/2023

% Plot smoothed fraction correct

% Plot each session running average of fraction correct

num_mice = size(tbt_cell,1);
num_days = size(tbt_cell,2);

figure
all_frac = nan.*ones(num_mice,num_days);
for m = 1:num_mice
    % cur_tbt = [];
    trial_num = zeros(num_days,1);
    for d = 1:num_days
        if ~isempty(tbt_cell{m,d})
            
            cur_tbt = tbt_cell{m,d}(3,ismember(tbt_cell{m,d}(3,:),7:12));
            all_frac(m,d) = sum(ismember(cur_tbt(1:av_num),[7,10]))/av_num;
            
        end        
    end
    
end

line_width = 2;
marker_size = 10;

xx = 1:num_days;
for m = 1:num_mice
    idx = ~isnan(all_frac(m,:));
    plot(xx(idx),squeeze(all_frac(m,idx)),'-o','LineWidth',line_width,'Color','k','MarkerFaceColor','k','MarkerSize',marker_size)
    hold on
end
box off
yticks([0,0.25,0.5,0.75,1])
ylim([0,1])
yline(0.5,'--','LineWidth',2);
xticks([1,2,3,4])
xlim([0.75,4.25])
title(["First 10 BMI Trial Performance"])
xlabel("Day")
ylabel("Fraction Correct")