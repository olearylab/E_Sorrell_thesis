function [xpred, Wout]  = Predict_and_update_ZA_LMS(zin, xin, W0, lrate, valid, update)
    persistent Win
    % ro = lrate*100;% original
    ro = 10^-3*lrate; % new
    if isempty(Win)
        Win = W0;
    end
    xpred = zin*Win;
    if valid && update
        error = xin - xpred;
        dW = lrate*zin'*error;
        % Wout = Win + dW - ro*sign(Win); % L1
        Wout = Win + dW - ro*Win; % Possibly L2?
        Win = Wout;
    else
        Wout = Win;
    end