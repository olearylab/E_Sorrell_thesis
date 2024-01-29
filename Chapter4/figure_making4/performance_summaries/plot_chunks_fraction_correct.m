function [] = plot_chunks_fraction_correct(tbt_cell,av_num)
% 10/03/2023

% Plot smoothed fraction correct

% Do I just plot a running average, or is he wanting me to smooth this?
% Smoothing just seems confounding and cheating?

num_mice = size(tbt_cell,1);
num_days = size(tbt_cell,2);

figure

for m = 1:num_mice
    cur_tbt = [];
    trial_num = zeros(num_days,1);
    for d = 1:num_days
        if ~isempty(tbt_cell{m,d})
            
            cur_tbt = [cur_tbt,tbt_cell{m,d}(3,ismember(tbt_cell{m,d}(3,:),7:12))];
            trial_num(d) = size(tbt_cell{m,d},2);
            
        end
    end
    
    running_res = zeros(ceil(length(cur_tbt)/av_num),1);
    for i = 1:ceil(length(cur_tbt)/av_num)
        if av_num*i<length(cur_tbt)
            cur_res = cur_tbt(av_num*(i-1)+1:av_num*i);
        else
            cur_res = cur_tbt(av_num*(i-1)+1:end);
        end
        running_res(i) = sum(ismember(cur_res,[7,10]))/length(cur_res);
    end
    
    plot(running_res,'LineWidth',2,'Color','k')
    hold on
    
end
ylim([0,1])
yline(0.5,'--','LineWidth',2);