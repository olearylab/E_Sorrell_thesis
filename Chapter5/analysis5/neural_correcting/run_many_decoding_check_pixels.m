function [all_res_cell,ball_weights_cell,bmi_weights_cell] = run_many_decoding_check_pixels(virmen_cell,train_mean_cell,Rfixed_cell,tbt_cell,model_params)
% 26/06/2023

%
% model_params.lrate = 10^-2;
% model_params.loops = 40;
% model_params.spatial = false;
% model_params.reg_images = false;
mice = ["DW81";"DW83";"DW86";"DW113";"DW129"];

data_days{1} = ["20200928";"20200929";"20200930";"20201001"];
data_days{2} = ["20200928";"20200929";"20200930";"20201001"];
data_days{3} = ["20200928";"20200929";"20200930";"20201001"];
data_days{4} = ["20210824";"NA";"20210826";"20210827"];
data_days{5} = ["20211117";"20211118";"20211119";"20211120"];

model_params.loops = 20;

num_mice = size(virmen_cell,1);
num_days = size(virmen_cell,2);

all_res_cell = cell(num_mice,num_days);
ball_weights_cell = cell(num_mice,num_days);
bmi_weights_cell = cell(num_mice,num_days);

for m = 1:num_mice
    data_d = data_days{m};
    for d = 1:num_days
        if ~isempty(virmen_cell{m,d})
            
            if m < 4
                res_path = '/Users/ethan/Dropbox (Cambridge University)/Boston_Data_09_2020/'; 
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
                cur_res_path = res_path + mice(m) + '/' + data_d(d) + '/test/';
                zfull = importdata(cur_res_path + mice(m) + '_' + data_d(d) + '_final_test.mat');
                zfull = permute_images(zfull);
            end
            disp("Running mouse " + m + " Day " + d)
            [all_res_mat,ball_weights,bmi_weights] = decoding_check(virmen_cell{m,d},zfull,tbt_cell{m,d},model_params,train_mean_cell{m},Rfixed_cell{m});
            all_res_cell{m,d} = all_res_mat;
            ball_weights_cell{m,d} = ball_weights;
            bmi_weights_cell{m,d} = bmi_weights;
            
        end
    end
end
