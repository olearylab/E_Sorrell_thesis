function [z_binned_res_cell,z_binned_cell,x_binned_cell,full_all_means] = get_z_binned_res_cell_samps_pix(virmen_cell,tbt_cell,Wout_online_cell,model_params,normalise_z,nbins,train_mean_cell,Rfixed_cell,use_bias)
% 22/09/2023
% For online view angle weights
% 

mice = ["DW81";"DW83";"DW86";"DW113";"DW129"];

data_days{1} = ["20200928";"20200929";"20200930";"20201001"];
data_days{2} = ["20200928";"20200929";"20200930";"20201001"];
data_days{3} = ["20200928";"20200929";"20200930";"20201001"];
data_days{4} = ["20210824";"NA";"20210826";"20210827"];
data_days{5} = ["20211117";"20211118";"20211119";"20211120"];

num_mice = size(virmen_cell,1);
num_days = size(virmen_cell,2);

z_binned_res_cell = cell(num_mice,num_days);
z_binned_cell = cell(num_mice,num_days);
x_binned_cell = cell(num_mice,num_days);

full_all_means = NaN(4,num_mice,num_days);

for m = 1:num_mice
    Wout = Wout_online_cell{m};
    Wout = Wout(:,2);
    data_d = data_days{m};
    for d = 1:num_days
        if ~isempty(virmen_cell{m,d})
            
            % x_vec = [6,7,15,17,on_ind,res_ind,vav_ind,vav_res_ind];
            % [z_binned] = calc_binned_residual_activity(z_cell{m,d},virmen_cell{m,d},tbt_cell{m,d},model_params,normalise_z,nbins);
            
            if m < 4
                res_path = '/Users/ethan/Dropbox (Cambridge University)/Boston_Data_09_2020/'; 
                % res_path = '/Users/ethansorrell/Dropbox (Cambridge University)/Boston_Data_09_2020/';
                cur_res_path = res_path + mice(m) + '/' + mice(m) + '_' + data_d(d) + '/test/';
                if isfile(cur_res_path + mice(m) + '_' + data_d(d) + '_final_test_revised.mat')
                    zfull = importdata(cur_res_path + mice(m) + '_' + data_d(d) + '_final_test_revised.mat');
                else
                    z1 = importdata(cur_res_path + mice(m) + '_' + data_d(d) + '_final_test_revised_1.mat');
                    z2 = importdata(cur_res_path + mice(m) + '_' + data_d(d) + '_final_test_revised_2.mat');
                    zfull = [z1;z2];
                end
            else
                res_path = '/Users/ethan/Dropbox (Cambridge University)/Boston_PPC_Data_05_2021/'; 
                % res_path = '/Users/ethansorrell/Dropbox (Cambridge University)/Boston_PPC_Data_05_2021/';
                cur_res_path = res_path + mice(m) + '/' + data_d(d) + '/test/';
                zfull = importdata(cur_res_path + mice(m) + '_' + data_d(d) + '_final_test.mat');
                zfull = permute_images(zfull);
            end
            disp("Running mouse " + m + " Day " + d)
            
            
            [z_binned,z_binned_res,x_binned,all_means] = calc_binned_residual_activity_from_samps_pix(zfull,virmen_cell{m,d},tbt_cell{m,d},Wout,model_params,normalise_z,nbins,train_mean_cell{m},Rfixed_cell{m},use_bias);

            z_binned_res_cell{m,d} = z_binned_res;
            z_binned_cell{m,d} = z_binned;
            x_binned_cell{m,d} = x_binned;
            full_all_means(:,m,d) = all_means;
                
        end
    end
end

save full_all_means_pixels_bias.mat full_all_means -v7.3
