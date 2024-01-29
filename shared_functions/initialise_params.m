% Inputs and parameters
% Use this script to set inputs and parameters. (Maybe not useful for
% iterating through different values?)
% Could save parameters with results for future reference.
% Could use separate variables as I have been in the past, or can use a
% struct.

%% struct model_params

% pre-processing
model_params.downsampled = true;
model_params.dff = true;
model_params.spatial = true;
model_params.reg_images = true;

% Learning rate - change to vector of rates for different variables?
if model_params.dff
    model_params.lrate = 10^-3;
else
    model_params.lrate = 10^-10;
end

% DFF offset
% 150 for original mice
% 200 for new mice 
model_params.dff_offset = 200;

% Filter constants
model_params.fs = 30;
model_params.blursmall = 0.6;
model_params.blurlarge = 5;
model_params.tau_fast = 0.5;
model_params.tau_slow = 45;
model_params.afast = (1-exp(-1/(model_params.tau_fast*model_params.fs))); 
model_params.aslow = (1-exp(-1/(model_params.fs*model_params.tau_slow)));

% Convergence criterion
model_params.loops = 5; % for fixed loops 

% Number of training trials
% model_params.num_correct_trials = 100; % not used

% Decoded variables
% [6,3,13,15,5,2,7,4] = [ypos,yvel,pitch,yaw,xpos,xvel,theta,thetavel]
model_params.xnums = [6,3,13,15,5,2,7,4];
% model_params.xnums = [13,15];

% Weight intialisation
model_params.Win = 0;

% Trial types for training
% Expand 1,2,3,4 into 1,4,7,10 
model_params.types_vec = [1,4];

% Optional z-scoring
model_params.zscore_z = false;
model_params.zscore_x = false;

% Optional registering to training data
model_params.train_reg = false;

% Optional view angle histogram balance
model_params.balance_hist = false;
model_params.balance_nbins = 6;

% If using offsets
model_params.subtract_offsets = false;
model_params.offsets = [1.4804,1.4844];

% If doing complex view angle decoding
model_params.va_2d = false;

% For balancing left and right in training
model_params.balance_lr = true;

% For regularising
model_params.regularise_training = false;
% If above is true, then true L1 is L1, false is L2
model_params.L1 = true;
% Regularisation parameter.
model_params.ro = 100*model_params.lrate;

% Optional removal of poor trials
model_params.remove_poor = true;
% Could also have parameter of view angle deemed poor.