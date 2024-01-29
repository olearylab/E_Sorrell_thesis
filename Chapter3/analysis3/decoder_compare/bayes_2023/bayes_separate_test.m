function [log_like_all,log_post_all, MAP_all] = bayes_separate_test(zfilt,xprev,all_priors,all_edges,all_means,varpix,nbins) 
% 26/05/2023

% Updated for many variables. 
% zimage should be filtered, ztestfilt

% Run Bayesian decoder on test data
num_x = size(all_means,3);

log_like_all = zeros(nbins,num_x);

for i = 1:nbins
    for j = 1:num_x
        log_like_all(i,j) = -sum(((all_means(:,i,j)'-zfilt).^2)./(2*varpix));
    end
end

% add prior using p(xt|xt-1)
prior_width_scale = 2; % 1 default
prior_strength_scale = 0.5; % 1 default
all_centres = (all_edges(:,2:end)+all_edges(:,1:end-1))/2;

MAP_all = zeros(num_x,1);
log_post_all = zeros(nbins,num_x);

for i = 1:num_x
    cur_pr = normpdf(all_centres(i,:),xprev(i),prior_width_scale*all_priors(i));

    log_post_all(:,i) = log_like_all(:,i) + log(cur_pr')*prior_strength_scale;
    
    [MAP,ind] = max(log_post_all(:,i));
    % MAP_all = (all_edges(ind)+all_edges(ind+1))/2;
    MAP_all(i) = all_centres(i,ind);

end


