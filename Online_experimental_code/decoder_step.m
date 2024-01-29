function [xpred]  = decoder_step(zin, W0)
%expecting to decode pitch and view angle
    persistent thetaprev
    dt = 1/30;
    beta = 4;
    xpred = zin*W0;
    if isempty(thetaprev) 
        thetaprev = xpred(2);
        xpred(2) = 0;
    else
        temp = xpred(2);
        xpred(2) = (xpred(2) - thetaprev)/(dt*-beta);
        thetaprev = temp;
    end
    end