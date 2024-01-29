function [d_image] = downsample_image(image, factor)
n = size(image,1)/factor;
d_image = zeros(n);
    for k = 1:n
        for l = 1:n
            cur_block = image(factor*(k-1)+1:factor*k,factor*(l-1)+1:factor*l);
            d_image(k,l) = squeeze(sum(cur_block,[1,2]))./(factor^2);
        end
    end
    d_image = d_image(:)';