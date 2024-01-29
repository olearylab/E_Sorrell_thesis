function [new_particles, pred, particle_weights, resampled_particles] = particle_general_sep_test(zfilt,particles,means,varpix,variances,edges)
% 26/05/2023

% Input data should already be filtered.

m = size(particles,1);
% Run Bayesian decoder on test data

log_like = zeros(size(particles));

%%%%% if binning particles %%%%%
binned = zeros(size(particles));
for i = 1:size(particles,2)
    binned(:,i) = discretize(particles(:,i),edges(i,:));
end

for i = 1:m
    for j = 1:size(particles,2)
        log_like(i,j) = 1/sum((((squeeze(means(j,:,binned(i,j)))-zfilt).^2)./(2*varpix)));
    end
%     log_like_y(i) = -sum(((meany(:,ybinned(i))'-zvar).^2)./(2*varpix)) - sum(log(sqrt(varpix))) - length(zvar)*log(sqrt(2*pi));
%     log_like_th(i) = -sum(((meanth(:,thbinned(i))'-zvar).^2)./(2*varpix)) - sum(log(sqrt(varpix))) - length(zvar)*log(sqrt(2*pi));
    %log_like(i) = prod(exp(-((meanyt(:,ybinned(i),thbinned(i))'-zvar).^2)./(2*varpix))./sqrt(2*pi*varpix));
end
%%%%%

%%%%% if interpolating %%%%%
% particle_means = zeros(size(means,1),size(means,2),m);
% for i = 1:size(particles,2)
%     centres = (edges(i,2:end)+edges(i,1:end-1))./2;
%     particle_means(i,:,:) = interp1(centres,squeeze(means(i,:,:))',particles(:,i))';
% end
% 
% for i = 1:m
%     for j = 1:size(particles,2)
%         log_like(i,j) = 1/sum((((squeeze(particle_means(j,:,i))-zvar).^2)./(2*varpix)));
%     end
% %     log_like_y(i) = -sum(((meany(:,ybinned(i))'-zvar).^2)./(2*varpix)) - sum(log(sqrt(varpix))) - length(zvar)*log(sqrt(2*pi));
% %     log_like_th(i) = -sum(((meanth(:,thbinned(i))'-zvar).^2)./(2*varpix)) - sum(log(sqrt(varpix))) - length(zvar)*log(sqrt(2*pi));
%     %log_like(i) = prod(exp(-((meanyt(:,ybinned(i),thbinned(i))'-zvar).^2)./(2*varpix))./sqrt(2*pi*varpix));
% end
%%%%%

% log_like_y(isnan(log_like_y)) = -inf;
% log_like_th(isnan(log_like_th)) = -inf;
log_like(isnan(log_like)) = 0;

%particle_weights_y = exp(log_like_y/length(zvar))/sum(exp(log_like_y/length(zvar)));
particle_weights = zeros(size(particles));
new_particles = zeros(size(particles));

for i = 1:size(particles,2)
    if (max(log_like(:,i))-min(log_like(:,i))) == 0
        particle_weights(:,i) = (log_like(:,i))/sum(log_like(:,i));
    else
        particle_weights(:,i) = (log_like(:,i)-min(log_like(:,i)))/sum(log_like(:,i)-min(log_like(:,i)));
    end
    particle_weights(isnan(particle_weights)) = 0;
    indices = datasample(1:m,m,'Weights',particle_weights(:,i));
    new_particles(:,i) = particles(indices,i);
end

pred = sum(new_particles)./m;
resampled_particles = new_particles;

var_scale = 1;
for i = 1:size(particles,2)
    %step_size = normrnd(pred(i)-prev_pred(i),sqrt(var_scale*variances(i)),[m,1]);
    step_size = normrnd(0,sqrt(var_scale*variances(i)),[m,1]);
    new_particles(:,i) = new_particles(:,i) + step_size;
    cur_particles = new_particles(:,i);
    cur_particles(cur_particles<edges(i,1)) = edges(i,1);
    cur_particles(cur_particles>edges(i,end)) = edges(i,end);
    new_particles(:,i) = cur_particles;
end
