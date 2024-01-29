%% set up parameters and constants
initialise_params;

%% Data handling
baseDir = '\\research.files.med.harvard.edu\Neurobio\HarveyLab\Tier1\Dan\BehaviorAnimals\Raw\';
animalName = 'DW129';
date = '20211117';
abfNumber = '002';
sessionNumber = '1';
abfFile = fullfile(baseDir,animalName,date,'abf',[date(3:end),abfNumber,'.abf']);
virmenFile = fullfile(baseDir,animalName,date,'virmen',date(3:end),['session_',sessionNumber],'sessionData.mat');


%%
zdat = permute(imgStack,[3 1 2]);
virmen_data = binnedVirmenData;
virmen_data = rearrange_virmen(virmen_data);

tic
disp('training start')
%[Wout, xpred, origztrainfilt, train_mean, Rfixed] = DW_online_training(zdat, virmen_data, xnums, lrate, loops, downsampled, dff, spatial, Win);
[Wout, xpred, origztrainfilt, train_mean, Rfixed, train_valid, model_params] = DW_online_training_05_2021(zdat, virmen_data, model_params);
save(fullfile(baseDir,animalName,date,'Wout.mat'),'Wout');
save(fullfile(baseDir,animalName,date,'model_params.mat'),'model_params');
for i =1:3;
    hSI.hUserFunctions.userFunctionsCfg(i).Enable = 1;
end
%
% Also save train_mean and Rfixed for image registration
save(fullfile(baseDir,animalName,date,'onlineImageData.mat'),'train_mean', 'Rfixed');

init_analogOut;


%
toc;
% disp('training complete')
