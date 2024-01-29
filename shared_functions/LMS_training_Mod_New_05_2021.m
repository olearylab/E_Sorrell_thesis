function [Wout, xpred, ztrainfilt, model_params] = LMS_training_Mod_New_05_2021(zimages, xvars, model_params, train_valid)
% Inner training function that pre-processes training data then performs
% LMS training

% zimages: samplesxheightxwidth (downsampled)
% if neurons, zimages: samplesxneurons
% xvars: samplesxvariables

% legacy parameter always true for LMS training
update = true;

train_length = size(xvars,1);

% allows for use with neurons or images
if ndims(zimages) == 3
    ztrainfilt = zeros(size(zimages,1),size(zimages,2)^2+1);
    % Initialise weight matrix
    if length(model_params.Win)==1
        W = zeros(size(zimages,2)^2+1,size(xvars,2));
        W(end,:) = mean(xvars,1);
    else
        W = model_params.Win;
    end
else
    ztrainfilt = zeros(size(zimages,1),size(zimages,2)+1);
    if length(model_params.Win)==1
        W = zeros(size(zimages,2)+1,size(xvars,2));
        W(end,:) = mean(xvars,1);
    else
        W = model_params.Win;
    end
end

% burn in dff filter (in reverse).
if model_params.dff
    for i = 1:size(xvars,1)
        % turn into array
        zvar = zimages((end+1)-i,:);

        % DF/F filtering
        [ztrainfilt(i,1:end-1)] = dff_filt(zvar, model_params.afast, model_params.aslow);

    end
end

% Now actually filter

for i = 1:size(xvars,1)
    % turn into array
    zvar = zimages(i,:);

    % DF/F filtering
    if model_params.dff
        [ztrainfilt(i,1:end-1)] = dff_filt(zvar, model_params.afast, model_params.aslow);
    else
        ztrainfilt(i,1:end-1) = zvar;
    end

    % Spatial bandpass filter of data
    if model_params.spatial
        [ztrainfilt(i,1:end-1)] = spatial_bandpass(ztrainfilt(i,1:end-1),model_params.blursmall,model_params.blurlarge);
    end
    ztrainfilt(i,end) = 1;
end

% room for editing here to accomodate for convergence criterion

% Optional z-scoring
if model_params.zscore_z
    zdim = size(ztrainfilt,2)-1;
    zmeans = zeros(zdim,1);
    z_vars = zeros(zdim,1);
    for n = 1:zdim
        zmeans(n) = mean(ztrainfilt(train_valid,n));
        z_vars(n) = var(ztrainfilt(train_valid,n));
        ztrainfilt(:,n) = (ztrainfilt(:,n) - zmeans(n))/sqrt(z_vars(n));
    end
    model_params.zmeans = zmeans;
    model_params.z_vars = z_vars;
end
    
if model_params.zscore_x
    xdim = size(xvars,2);
    xmeans = zeros(xdim,1);
    x_vars = zeros(xdim,1);
    for n = 1:xdim
        xmeans(n) = mean(xvars(train_valid,n));
        x_vars(n) = var(xvars(train_valid,n));
        xvars(:,n) = (xvars(:,n) - xmeans(n))/sqrt(x_vars(n));
    end
    model_params.xmeans = xmeans;
    model_params.x_vars = x_vars;
    W(end,:) = mean(xvars,1);
end

% Optional offset subtraction
if model_params.subtract_offsets && sum(ismember(model_params.xnums,13))>0
    xvars(:,model_params.xnums==13) = xvars(:,model_params.xnums==13) - model_params.offsets(1);
    if sum(ismember(model_params.xnums,15))>0 && ~model_params.balance_hist
        xvars(:,model_params.xnums==15) = xvars(:,model_params.xnums==15) - model_params.offsets(2);
    end
end

for i = 1:model_params.loops
    % randomly shuffle training data
    shuffle = randperm(train_length);

    for j = 1:train_length
        
        if model_params.regularise_training
            % Can replace with training using regularisation e.g. ZA (L1)
            [xpred, Wout]  = Predict_and_update_ZA_LMS_2022(ztrainfilt(shuffle(j),:), xvars(shuffle(j),:), W, model_params, train_valid(shuffle(j)), update);
        else
            % LMS
            [xpred, Wout]  = Predict_and_update(ztrainfilt(shuffle(j),:), xvars(shuffle(j),:), W, model_params.lrate, train_valid(shuffle(j)), update);
        end
        % Optional plotting of weights during training
        % im = reshape(Wout(1:end-1,1),128,128);
        % imagesc(im)
        % drawnow
    end
    % Optional plotting of weight image after each loop
%     subplot(ceil(model_params.loops/ceil(sqrt(model_params.loops))),ceil(sqrt(model_params.loops)),i)
%     im = reshape(Wout(1:end-1,1),128,128);
%     imagesc(im)
%     title("Loop " + i)
%     drawnow
end

% Optional va histogram balancing
% could do before all training, but this would effect samples for pitch
% training as well (which we don't want, or may want to balance
% separately).
clear Predict_and_update
if model_params.balance_hist
    if sum(ismember(model_params.xnums,7))>0
        [wrapped_va_new,ztrainfilt_new] = balance_training_histogram_general_flat(xvars(train_valid,model_params.xnums==7),ztrainfilt(train_valid,:),model_params.balance_nbins);
        train_length = length(wrapped_va_new);
        % only used valid data, so all data now valid.
        train_valid_new = ones(train_length,1);
        % set mean back to 0
        W(end,model_params.xnums==7) = 0;
        for i = 1:model_params.loops
        % randomly shuffle training data
        shuffle = randperm(train_length);

            for j = 1:train_length
                % LMS
                [xpred, Wout_new]  = Predict_and_update(ztrainfilt_new(shuffle(j),:), wrapped_va_new(shuffle(j)), W(:,model_params.xnums==7), model_params.lrate, train_valid_new(shuffle(j)), update);

                % Can replace with training using regularisation e.g. ZA (L1)
                % [xpred, Wout]  = Predict_and_update_ZA_LMS(ztrainfilt(shuffle(j),:), xvars(shuffle(j),:), W, model_params.lrate, train_valid(shuffle(j)), update);
            end
        end
        Wout(:,model_params.xnums==7) = Wout_new;
    end
    clear Predict_and_update
    if sum(ismember(model_params.xnums,15))>0
        offset = model_params.offsets(2);
        [yaw_new,ztrainfilt_new] = balance_training_histogram_general_flat(xvars(train_valid,model_params.xnums==15)-offset,ztrainfilt(train_valid,:),model_params.balance_nbins);
        if model_params.subtract_offsets == false
            yaw_new = yaw_new + offset;
        end
        train_length = length(yaw_new);
        % only used valid data, so all data now valid.
        train_valid_new = ones(train_length,1);
        % set mean to new mean
        W(end,model_params.xnums==15) = mean(yaw_new);
        for i = 1:model_params.loops
        % randomly shuffle training data
        shuffle = randperm(train_length);

            for j = 1:train_length
                % LMS
                [xpred, Wout_new]  = Predict_and_update(ztrainfilt_new(shuffle(j),:), yaw_new(shuffle(j)), W(:,model_params.xnums==15), model_params.lrate, train_valid_new(shuffle(j)), update);

                % Can replace with training using regularisation e.g. ZA (L1)
                % [xpred, Wout]  = Predict_and_update_ZA_LMS(ztrainfilt(shuffle(j),:), xvars(shuffle(j),:), W, model_params.lrate, train_valid(shuffle(j)), update);
            end
        end
        Wout(:,model_params.xnums==15) = Wout_new;
    end
end