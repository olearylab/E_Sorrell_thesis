function [pix_neuron_list,pix_neuron_mat,sorted_pix_list,sorted_pix_is_neuron,final_footprints] = assign_pixel_weights_to_neurons_11042023(Wout,stat)
% 28/09/2022

% Attempt to assign each pixel any neurons that have masks that overlap
% with 
% Wout is input as 128x128 in correct transpose (i.e. dealt with outside)

% Calculate neurons weights from pixels using neuron masks.
% Calculate both max and mean. Can decide which to use later
% Need to downsample neruon masks... how does this work??? what method to
% use?
% Need to be careful about dimensions, need to make sure mask matches
% weights. Might be different old vs new mice? Not sure

num_neurons = length(stat);
% num_x = size(Wout,2);
% n_weights_max_abs = zeros(num_neurons,num_x);
% n_weights_mean = zeros(num_neurons,num_x);
%% 
pix_neuron_mat = cell(128,128);
pix_neuron_list = cell(128*128,1);

final_footprints = zeros(128,128);
for n = 1:num_neurons
    cur_stat = stat{n};
    cur_footprint = zeros(512,512);
    for j = 1:length(cur_stat.xpix)
        cur_footprint(cur_stat.xpix(j)+1,cur_stat.ypix(j)+1) = 1;
    end
    down_footprint = zeros(128,128);
    for i = 1:128
        for j = 1:128
            if sum(cur_footprint((i-1)*4+1:i*4,(j-1)*4+1:j*4),'all')>0
                down_footprint(i,j) = 1;
                final_footprints(i,j) = 1;
            end
        end
    end
    
    for i = 1:128
        for j = 1:128
            if down_footprint(i,j) == 1
                pix_neuron_mat{i,j} = [pix_neuron_mat{i,j},n];
            end
        end
    end
    
    down_foot_flat = down_footprint(:);
    for i = 1:length(down_foot_flat)
        if down_foot_flat(i) == 1
            pix_neuron_list{i} = [pix_neuron_list{i},n];
        end
    end
   
%     
%     for i = 1:num_x
%         cur_w = reshape(Wout(1:end-1,i),128,128);
%         cur_n_w = cur_w(down_footprint==1);
%         if ~isempty(cur_n_w)
%             n_weights_max_abs(n,i) = max(abs(cur_n_w));
%             n_weights_mean(n,i) = mean(cur_n_w);
%         end
%     end
    
end

for i = 1:128
    for j = 1:128
        if isempty(pix_neuron_mat{i,j})
            pix_neuron_mat{i,j} = 0;
        end
    end
end

for i = 1:(128*128)
    if isempty(pix_neuron_list{i})
        pix_neuron_list{i} = 0;
    end
end
%% Do something with weights
% Am I certain weights match masks (is there a transpose needed?) This
% might be different between those trained online vs trained offline?)

% Wout is 128x128

% Output sorted list of pixels by weights largest to smallest absolute
% values
sorted_pix_list = cell(128*128,1);
[sorted_w,w_inds] = sort(abs(Wout(:)),'descend');
sorted_pix_is_neuron = zeros(128*128,1);
for i = 1:length(pix_neuron_list)
    sorted_pix_list{i} = pix_neuron_list{(w_inds(i))};
    sorted_pix_is_neuron(i) = sum(sorted_pix_list{i}) > 0;
end