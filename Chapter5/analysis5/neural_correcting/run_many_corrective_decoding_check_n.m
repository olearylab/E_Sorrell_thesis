function [all_res_cell,all_samps] = run_many_corrective_decoding_check_n(virmen_cell,z_cell,tbt_cell,virmen_train_cell,tbt_train_cell,model_params,offsets,num_subs)
% 08/09/2023
nbins = 50;
%
model_params.lrate = 10^-2;
model_params.loops = 20;
model_params.spatial = false;
model_params.reg_images = false;

num_mice = size(virmen_cell,1);
num_days = size(virmen_cell,2);

all_res_cell = cell(num_mice,num_days);
all_samps = NaN(num_mice,num_days,3);

for m = 1:num_mice
    [mean_binned,std_binned] = calculate_mean_ball_va(virmen_train_cell{m},tbt_train_cell{m},nbins,offsets(m,:));
    for d = 1:num_days
        if ~isempty(z_cell{m,d})
            disp("Running mouse " + m + " Day " + d)
            % [all_res_mat] = decoding_check(virmen_cell{m,d},z_cell{m,d},tbt_cell{m,d},model_params,[],[]);
            [all_res_mat,bmi_weights,num_samps] = decoding_check_corrective(virmen_cell{m,d},z_cell{m,d},tbt_cell{m,d},model_params,[],[],nbins,offsets(m,:),mean_binned,std_binned,num_subs);
            all_res_cell{m,d} = all_res_mat;
            all_samps(m,d,:) = num_samps;
        end
    end
end
