function [accuracies_mat] = svm_behaviour_classifier_15052023(s)
% 15/05/2023

% Updated SVM classification of behaviour.

% I often have very few trials to train and test with. Since this is
% behaviour not neural activity, maybe I should try combining multiple
% sessions?

% Maybe try to make it all as similar as possible to the LSTM (Same
% training testing split etc. This requires doing 5 fold cross validation
% then).

% Should I do bin specific or just general? - bin specific

% edit for when decoding vr pitch and view angle
p_ind = size(s.binnedVirmenData2,2);
% tempVirmenData = s.binnedVirmenData2(:,[p_ind,7],:);

tempVirmenData = s.binnedVirmenData2(:,[13,15],:);
% tempVirmenData = s.binnedVirmenData2(:,[13,15,7],:);
tempVirmenData(:,1,:) = zscore(tempVirmenData(:,1,:),0,'all');
tempVirmenData(:,2,:) = zscore(tempVirmenData(:,2,:),0,'all');
% tempVirmenData(:,3,:) = zscore(tempVirmenData(:,3,:),0,'all');

nbins = size(tempVirmenData,3);

c_type = 'svm';
% Don't trust inbuilt cross-validation, from some other classifying I did.

% for trialNumber =1:size(s.binnedVirmenData2,1);
%     binnedTrialData{trialNumber} = squeeze(tempVirmenData(trialNumber,:,:));
% end

% set up partitions
% I think this is how to train only on correct trials
s.controlIdx = find(s.trialType<3);
if s.correct_only
    s.controlIdx = s.controlIdx(s.correctVec(s.controlIdx)==1);
end
s.controlBinnedData = tempVirmenData(s.controlIdx,:,:);
s.controlLabels = s.trialType(s.controlIdx);

s.controlLabelSequence = zeros(size(s.controlBinnedData,1),size(s.controlBinnedData,3));

for controlTrialNumber = 1:size(s.controlBinnedData,1);
    s.controlLabelSequence(controlTrialNumber,:) = categorical(repmat(s.controlLabels(controlTrialNumber),[1, size(tempVirmenData,3)]),[1 2]);
end

nModels = 5; % number of cv partitions

c = cvpartition(size(s.controlLabelSequence,1)/2,'KFold',nModels);
testSplitID = NaN(size(s.trialType));

% Save training and validation accuracies
% accuracies_mat = zeros(2,nbins);

many_Mdls = cell(nModels,nbins);
% predict_vec = zeros(size(s.controlBinnedData,1),nbins);

% Trying general classification first, then bin specific.
% Changed to bin specific now.

%% Train network
for cvSplitNumber=1:nModels
trainIdx = [c.training(cvSplitNumber);c.training(cvSplitNumber)];
testIdx = [c.test(cvSplitNumber);c.test(cvSplitNumber)];
   
XTrain = s.controlBinnedData(trainIdx,:,:);
YTrain = s.controlLabelSequence(trainIdx,:);

XTest = s.controlBinnedData(testIdx,:,:);
YTest = s.controlLabelSequence(testIdx,:);

XTrain = permute(XTrain,[1,3,2]);
XTest = permute(XTest,[1,3,2]);
% XTrain = XTrain(:,:)';
% YTrain = YTrain(:);

tidx = find(testIdx==1);
for i = 1:nbins
    many_Mdls{cvSplitNumber,i} = fitclinear(squeeze(XTrain(:,i,:)),YTrain(:,i),'Learner',c_type);
    % many_Mdls{cvSplitNumber,i} = fitclinear(squeeze(XTrain(:,i,:)),YTrain(:,i),'OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',struct('Verbose',0,'ShowPlots',false));
    
    % many_Mdls{cvSplitNumber,i} = fitclinear(squeeze(XTrain(:,i,:)),YTrain(:,i),'Learner',c_type,'Lambda',lambda);
    % bin_test = predict(many_Mdls{cvSplitNumber,i},squeeze(XTest(:,i,:)));
    % predict_vec(testIdx,i) = bin_test;
end

% testSplitID(s.controlIdx(c.test(cvSplitNumber)))=1;
testSplitID(s.controlIdx(testIdx))=cvSplitNumber;
end

predict_vec = zeros(length(s.trialType),nbins);
accuracies_mat = zeros(4,nbins);
correct_mat = zeros(length(s.trialType),nbins);

testTrialTypes = s.trialType;
testTrialTypes(testTrialTypes>2) = testTrialTypes(testTrialTypes>2) - 2;

for trialNumber= 1:length(s.trialType);
if ~isnan(testSplitID(trialNumber));
    for i = 1:nbins
        bin_test = predict(many_Mdls{testSplitID(trialNumber),i},squeeze(tempVirmenData(trialNumber,:,i)));
        predict_vec(trialNumber,i) = bin_test;
        correct_mat(trialNumber,i) = predict_vec(trialNumber,i)==testTrialTypes(trialNumber);
    end
else
    temp1 = [];
    for i = 1:nbins
        for modelNumber = 1:nModels
            bin_test = predict(many_Mdls{modelNumber,i},squeeze(tempVirmenData(trialNumber,:,i)));
            temp1 = [temp1;bin_test==testTrialTypes(trialNumber)];
        end
        correct_mat(trialNumber,i) = squeeze(nanmean(temp1,1));
    end
end
end;

% correct_mat = predict_vec==s.controlLabelSequence;
% accuracies_mat(1,:) = sum(correct_mat(s.controlLabels==1,:),1)./sum(s.controlLabels==1);
% accuracies_mat(2,:) = sum(correct_mat(s.controlLabels==2,:),1)./sum(s.controlLabels==2);

for i = 1:4
    accuracies_mat(i,:) = sum(correct_mat(s.trialType==i & s.correctVec==1,:),1)./sum(s.trialType==i & s.correctVec==1);
end

% plot(accuracies_mat')