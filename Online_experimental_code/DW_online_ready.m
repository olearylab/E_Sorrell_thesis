%% flags
downsampled = true;
dff = true;
spatial = true;
update = true;
loops = 5;
% Need to load in a saveed Wout and assign to Win
load('\\research.files.med.harvard.edu\Neurobio\HarveyLab\Tier1\Dan\BehaviorAnimals\Raw\DW83\20200925\Wout.mat');
Win = Wout;

%xnums = [6,3,13,15,5,2,7,4];
xnums = [13,7];
lrate = 10^-10;
if dff
    lrate = 10^-3;
end
fs = 30;

afast = (1-exp(-1/(0.5*fs)));
aslow = (1-exp(-1/(fs*45)));
%lrate = 10^-10;
blursmall = 0.6;
blurlarge = 5;

%% Data handling
baseDir = '\\research.files.med.harvard.edu\Neurobio\HarveyLab\Tier1\Dan\BehaviorAnimals\Raw\';
animalName = 'DW83';
date = '20201002';
abfNumber = '004';
sessionNumber = '1';
abfFile = fullfile(baseDir,animalName,date,'abf',[date(3:end),abfNumber,'.abf']);
virmenFile = fullfile(baseDir,animalName,date,'virmen',date(3:end),['session_',sessionNumber],'sessionData.mat');


%%
zdat = permute(imgStack,[3 1 2]);
virmen_data = binnedVirmenData;

tic
disp('training start')
[Wout, xpred, origztrainfilt, train_mean, Rfixed] = DW_online_training(zdat, virmen_data, xnums, lrate, loops, downsampled, dff, spatial, Win);
save(fullfile(baseDir,animalName,date,'Wout.mat'),'Wout');
for i =1:3;
    hSI.hUserFunctions.userFunctionsCfg(i).Enable = 1;
end
%
save(fullfile(baseDir,animalName,date,'onlineImageData.mat'),'train_mean', 'Rfixed');

init_analogOut;
% Also save train_mean and Rfixed for image registration

%
toc;
% disp('training complete')
% xprediction = zeros(size(xtest,1),length(xnums));



% tic;
% disp('testing start')
% for i = 1:size(xtest,1)
% [xpred] = Decoder_Online_Mod(squeeze(ztest(i,:,:)), Wout, downsampled, afast, aslow, blurlarge, blursmall, dff, spatial);
% xprediction(i,:) = xpred;
% end
% toc;
% disp('testing complete')