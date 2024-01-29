function [means,varpix,edges,variances] = particle_general_sep_train(ztrainfilt,xtrain,nbins,xnums)
% 26/05/2023

% Updated. Input is filtered zdata

% Estimate distributions for particle filter
n = size(ztrainfilt,2); % now 16348

% get values and remove ITI
variables = xtrain(xnums,:);

ITI = xtrain(8,:);
[train_valid] = clean_valid_data(ITI);
variables = variables(:,train_valid);

trials = xtrain(12,:);
trials = trials(train_valid);
varsteps = [];
varacc = [];
for i = 1:trials(end)
    curvars = variables(:,trials==i);
    varsteps = [varsteps, curvars(:,2:end)-curvars(:,1:end-1)];
    %for momentum prior
    varacc = [varacc,varsteps(:,2:end)-varsteps(:,1:end-1)];
end

variables(xnums==7,:) = wrapToPi(variables(xnums==7,:));

% [N,Xedges,Yedges] = histcounts2(y,th,nbins);
edges = zeros(length(xnums),nbins+1);
binned = zeros(length(xnums),size(variables,2));
for i = 1:length(xnums)
    [~,edges(i,:)] = histcounts(variables(i,:),nbins);
    binned(i,:) = discretize(variables(i,:),nbins);
end

ztrainfilt = ztrainfilt(train_valid,:);

means = zeros(length(xnums),n,nbins);

for i = 1:nbins
    for j = 1:length(xnums)
        means(j,:,i) = mean(ztrainfilt(binned(j,:)==i,:));
    end
end
varpix = var(ztrainfilt);

%covyt = [mean(abs(ysteps))^2,0;0,mean(abs(thsteps))^2];

% for velocity prior
variances = max(abs(varsteps)').^2;
%variances = mean(abs(varsteps)').^2;
% for momentum prior
%variances = max(abs(varacc)').^2;