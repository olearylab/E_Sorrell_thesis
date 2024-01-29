function [xpred]  = decoder_step(zin, W0)
    xpred = zin*W0;
    end