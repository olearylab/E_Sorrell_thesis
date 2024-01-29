function [all_means,varpix,all_priors,all_edges] = bayes_separate_train(ztrainfilt,xtrain,nbins,xnums)
% 26/05/2023

% Updated to allow for more variables.
% Assumes training data is already preprocessed.

% Train Bayesian decoder, i.e. calculate prior and relevant means and
% variances.
n = size(ztrainfilt,2);

% xtrain(7,:) = wrapToPi(xtrain(7,:));

% get values and remove ITI
xvars = xtrain(xnums,:);
ITI = xtrain(8,:);
[train_valid] = clean_valid_data(ITI);

xvars = xvars(:,train_valid);

% alternate prior using p(xt|xt-1) as gaussian
trials = xtrain(12,:);
trials = trials(train_valid);
all_steps = [];

for i = 1:trials(end)
    cur_vars = xvars(:,trials==i);

    all_steps = [all_steps, cur_vars(:,2:end)-cur_vars(:,1:end-1)];
end
all_priors = mean(abs(all_steps),2,'omitnan');

% wrap after calculation of step sizes.
xvars(xnums==7,:) = wrapToPi(xvars(xnums==7,:));

% calculate edges
all_edges = nan.*ones(length(xnums),nbins+1);
for i = 1:length(xnums)
    [~,all_edges(i,:)] = histcounts(xvars(i,:),nbins);
end

% bin data
all_binned = nan.*ones(size(xvars));
for i = 1:length(xnums)
    all_binned(i,:) = discretize(xvars(i,:),nbins);
end


ztrainfilt = ztrainfilt(train_valid,:);

all_means = nan.*ones(n,nbins,length(xnums));

for i = 1:nbins
    for j = 1:length(xnums)
        all_means(:,i,j) = mean(ztrainfilt(all_binned(j,:)==i,:));
    end
end

% Maybe should be bin specific
varpix = var(ztrainfilt);

all_vars = nan.*ones(n,nbins,length(xnums));

for i = 1:nbins
    for j = 1:length(xnums)
        all_vars(:,i,j) = var(ztrainfilt(all_binned(j,:)==i,:));
    end
end