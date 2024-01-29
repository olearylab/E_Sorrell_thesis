function [] = plot_smoothed_fraction_correct_thesis(tbt_cell,av_num)
% 09/06/2023

% Plot smoothed fraction correct

% Plot each session running average of fraction correct

num_mice = size(tbt_cell,1);
num_days = size(tbt_cell,2);

% num_days = 1;

figure

for m = 1:num_mice
    % cur_tbt = [];
    trial_num = zeros(num_days,1);
    for d = 1:num_days
        if ~isempty(tbt_cell{m,d})
            
            cur_tbt = tbt_cell{m,d}(3,ismember(tbt_cell{m,d}(3,:),7:12));
            trial_num(d) = size(tbt_cell{m,d},2);          
            
            running_res = zeros(length(cur_tbt)-(av_num-1),1);
            for i = 1:length(cur_tbt)-(av_num-1)
                cur_res = cur_tbt(i:i+av_num-1);
                running_res(i) = sum(ismember(cur_res,[7,10]))/length(cur_res);
            end

            plot(running_res,'k','LineWidth',2)
            % plot(movmean(ismember(cur_tbt,[7,10]),av_num),'k','LineWidth',2)
            hold on
        end        
    end
    
end

ylim([0,1])
yline(0.5,'--','LineWidth',2);