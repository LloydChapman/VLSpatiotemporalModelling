%RUNPOSTPROCESSING Script to process MCMC output

%% Run script containing parameter values for different models
SpecifyModel

%% Set parameters for processing output
burnin=round(niters/5); % burn-in to discard
str='DIC'; % string to prepend to rslts when saving DIC
ss='mode'; % summary statistic to use in calculation of DIC
% str='DIC_mean'; % string to prepend to rslts when saving DIC
% ss='mean'; % summary statistic to use in calculation of DIC
np=7; % number of parameters that can be estimated

%% Calculate DIC
RunCalcDIC(rslts,burnin,ss,str)

%% Output results and find best-fitting model (with lowest DIC)
DICminMdl=RunCalcParEstsAndDICdiffs2(rslts,np,burnin,IPD,str,true,false,true,runName);

%% Plot deviance distributions
ord=[3,1,5,4,2,6];
RunPlotDevDistn(rslts(ord),burnin,IPD,h0s(ord),h40s(ord),delta0s(ord))

%% Calculate transmission probabilities for best-fitting model
thin=1;
CalcTrnsmssnProbs(rslts{DICminMdl},burnin,thin,runName)

%% Calculate contributions of different infectious groups for best-fitting model
nsmpls=min(1000,niters-burnin);
[Mcntrbtn,HPDIcntrbtn,iters]=CalcCntrbtnMgrtnAsx(rslts{DICminMdl},nsmpls,burnin);

%% Calculate reproduction numbers and infection distances and intervals for best-fitting model
[infctn,infctr,src,dists,times,SI,onsetinfctr,rcvryinfctr,meandists,meantimes,infctrmax,srcmax,distsmax,timesmax,SImax,onsetinfctrmax,rcvryinfctrmax,meandistsmax,meantimesmax,RjA_I,RjI_I,RjP_I,Rj_I,Rts_I,Rt_I,diA_I,diI_I,diP_I,di_I,tiA_I,tiI_I,tiP_I,ti_I,infctnA,infctrA,srcA,distsA,timesA,SIA,onsetinfctrA,rcvryinfctrA,infctrmaxA,srcmaxA,distsmaxA,timesmaxA,SImaxA,onsetinfctrmaxA,rcvryinfctrmaxA,RjA_A,RjI_A,RjP_A,Rj_A,Rts_A,Rt_A,diA_A,diI_A,diP_A,di_A,tiA_A,tiI_A,tiP_A,ti_A,iters,OT,Rj,Rts,Rt]=CalcInfctrMgrtnAsx(rslts{DICminMdl},nsmpls,burnin,iters);

%% Plot consensus transmission tree for part of the study area in the SE cluster for the best-fitting model
lon_rng=[90.295 90.303];
lat_rng=[24.609 24.617];
PlotTrnsmssnTreeCnsnssAsx(rslts{DICminMdl},nsmpls,infctn,infctr,src,infctnA,infctrA,srcA,iters,burnin,lon_rng,lat_rng)

%% Plot transmission tree from single sample for same part of the study area
[psmpl,p1smpl]=PlotTrnsmssnTreeSnglSmplAsx(rslts{DICminMdl},infctn,infctr,src,infctnA,infctrA,srcA,iters,burnin,lon_rng,lat_rng);

%% Reconstruct transmission tree for best-fitting model
[gen,dist,time,chain]=ReconstructTrnsmssnTreeCnsnss(rslts{DICminMdl},[infctn;infctnA],[infctr;infctrA],[src;srcA],infctnA);
save('TrnsmssnTree','gen','dist','time','chain')