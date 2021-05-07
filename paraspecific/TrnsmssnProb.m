function [Mprob,HPDIprob]=TrnsmssnProb(px,p,hx,nbins)
prob=1-exp(-px.*p*hx);
[Mprob,HPDIprob]=CalcModeAndHPDI(prob,nbins);