function [xvar_new,ztrain_new] = balance_training_histogram_general(xvar,ztrain,nbins)

% balance positive and negaitve histogram samples for training to
% reduce/prevent decoder bias.

% need to determine number of bins we care about? anything from 10-100
% probably ok?

% input is any training variable. If yaw, offset must already be
% subtracted, then add back on after... This is dumb, could just keep
% subtracted then not subtract in testing. This is much smarter!

% bin data/ create histogram
% create bins encompassing full range, and symmetric about 0
% xvar = xvar;
max_val = max(abs(xvar));
edges = linspace(-max_val,max_val,nbins+1);
[binned] = discretize(xvar,edges);

% should be equal number of bins either side of zero. Assume even nbins.
% Issue of 0 in some bins on one side? Either have to ignore, or change
% bins, or remove all samples from corresponding bin from other direction.
% For now just ignore but print warning
indices = [];
for i = 1:(nbins/2)
    sum_bin_neg = sum(binned == i);
    sum_bin_pos = sum(binned == (nbins - (i-1)));
    if (sum_bin_neg == 0)
        disp("Warning: Bin " + i + " has 0 samples") 
    elseif (sum_bin_pos == 0)
        disp("Warning: Bin " + (nbins - (i-1)) + " has 0 samples") 
    elseif sum_bin_neg > sum_bin_pos
        num_extra = sum_bin_neg - sum_bin_pos;
        cur_binned = find(binned==(nbins - (i-1)));
        for j = 1:floor(num_extra/sum_bin_pos)
            % wrapped_va_new = [wrapped_va_new;wrapped_va(cur_binned)];
            % ztrain_new = [ztrain_new;ztrain(cur_binned,:)];
            indices = [indices;cur_binned];
        end
        extras_vec = randperm(sum_bin_pos,rem(num_extra,sum_bin_pos));
        % wrapped_va_new = [wrapped_va_new;wrapped_va(cur_binned(extras_vec))];
        % ztrain_new = [ztrain_new;ztrain(cur_binned(extras_vec),:)];
        indices = [indices;cur_binned(extras_vec)];
    elseif sum_bin_neg < sum_bin_pos
        num_extra = sum_bin_pos - sum_bin_neg;
        cur_binned = find(binned==i);
        for j = 1:floor(num_extra/sum_bin_neg)
            % wrapped_va_new = [wrapped_va_new;wrapped_va(cur_binned)];
            % ztrain_new = [ztrain_new;ztrain(cur_binned,:)];
            indices = [indices;cur_binned];
        end
        extras_vec = randperm(sum_bin_neg,rem(num_extra,sum_bin_neg));
        % wrapped_va_new = [wrapped_va_new;wrapped_va(cur_binned(extras_vec))];
        % ztrain_new = [ztrain_new;ztrain(cur_binned(extras_vec),:)];
        indices = [indices;cur_binned(extras_vec)];
    end
end

xvar_new = [xvar;xvar(indices)];
ztrain_new = [ztrain;ztrain(indices,:)];

% figure(1)
% histogram(xvar,edges)
% figure(2)
% histogram(xvar_new,edges)
