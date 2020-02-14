function RunMCMC(i)
%RUNMCMC Code to select model and run MCMC
%   RUNMCMC(I) runs MCMC for model I as specified in SPECIFYMODEL and for 
%   niters iterations as set in SPECIFYMODEL.

%% Run script which loads data and sets parameter values for different models
SpecifyModel

% Set maximum number of computational threads according to number of 
% threads available (here 16) and number of models run (here 6), e.g. as
% floor(#threads/#models)
maxNumCompThreads(2);

%% Run MCMC
VLStmMCMC2010FastLklhdMgrtn(data,r1,p10,a,b,p2,u{i},beta0s(i),alpha0s(i),epsilon0,delta0s(i),lambda0s(i),h0s(i),h1s(i),h2s(i),h3s(i),hmssngs(i),h40s(i),pI0s(i),typ{i},inclLST(i),niters,plotOutpt,rslts{i},para)