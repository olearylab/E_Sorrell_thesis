function [xpred]  = decoder_step_basic(zin, W0)
    xpred = zin*W0;
    end