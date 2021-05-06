%RUNDICTEST Script to test DIC on simulated data
clear

%% DATA
% Load data
load('data_final2.mat'); % loads table called 'data'

% Set which paras are included in data
para=1:4; %1
str='_Paras1to4';

%% MODEL PARAMETERS
% Parameters that can be estimated
% 1: beta - spatial transmission rate constant
% 2: alpha - distance scale factor for spatial kernel
% 3: epsilon - background transmission rate
% 4: delta - additional within-household transmission rate
% 5: lambda0 - historical asymptomatic infection rate
% 6: h4 - relative infectiousness of asymptomatics
% 7: pI - proportion of infections that lead to VL

% Set parameters for negative binomial NB(r1,p1) incubation period distribution
r1=3; % shape parameter
mu=5; % guess for mean incubation period
p10=r1/(mu-1+r1); % initial guess for 'success' parameter
% Shape parameters for beta prior for p1
a=22;
b=(mu-1)*(a-1)/r1; % since mu-1=r1*b/(a-1)

% Set parameter for geometric asymptomatic infection period distribution
p2=1/5;

% Set parameters to estimate
u=1:4;

% Initial guesses/fixed values for transmission parameters
beta0=3;
alpha0=100;
epsilon0=1e-3;
delta0s=[0,1e-2];
ndelta0s=numel(delta0s);
pI0=0.15;
% Estimate lambda0 by fitting initial status model to LST data from 2002
s1=load('data_final');
[pars,~]=FitCatModLST3(s1.data,p2);
lambda0=pars/12;

% Relative infectiousnesses of different forms of PKDL
h1=9/26/(10/15); % macular/papular
h3=18/21/(10/15); % nodular
h2=(h1+h3)/2; % plaque (assumed halfway between macular/papular and nodular)
hmssng=(101/138*h1+31/138*h2+6/138*h3); % unexamined - use average relative infectiousness of examined PKDL cases for cases who weren't physically examined
h40=0.02; % relative infectiousness of asymptomatics and pre-symptomatics

% Set kernel type
typ='Exp';

%% MCMC PARAMETERS
niters=1e5; % number of iterations
plotOutpt=false; % flag for whether to plot output in real-time
if r1>1
    IPD='NBIP';
elseif r1==1
    IPD='GIP';
end

%% RUN MCMC ON SIMULATIONS WITHOUT AND WITH ADDITIONAL WITHIN-HOUSEHOLD TRANSMISSION
runName=[str '_DIC_test_addtnl_wthnHH_trnsmssn'];
% Names of models for saving
Mdls={'_HGHR_PRESX_ASX_INFCTS','_DELTA_HGHR_PRESX_ASX_INFCTS'};
rslts=cellfun(@(x)['MCMC_' IPD '_PKDL_ASX' x runName],Mdls,'UniformOutput',false);

% Run MCMC for each model on each of the simulated datasets without and
% with additional within-household transmission
% Make matrix of all combinations of indices for additional within-household 
% transmission in simulations and in models
[X,Y]=meshgrid(1:ndelta0s,1:ndelta0s);
Z=[X(:),Y(:)];
delta0sc=delta0s(Z(:,2));
uc=cell(1,ndelta0s^2);
uc(delta0sc==0)={1:3};
uc(delta0sc~=0)={1:4};
rsltsc=strcat(rslts(Z(:,2)),"_",cellstr(num2str(Z(:,1)))');
parfor i=1:ndelta0s^2
    % Set maximum number of computational threads for each parallel worker
    maxNumCompThreads(floor(12/ndelta0s^2))
    % Run each model on each of the simulated datasets (j) without and with additional within-household transmission
    j=Z(i,1);
    % Load simulated data
    simdata=readtable(['sim_output' str '_addtnl_wthnHH_trnsmssn_' num2str(j) '.csv'],'ReadVariableNames',true); % simulated data for para 1 with ith presymptomatic/asymptomatic infectiousness
    Status0=csvread(['init_status' str '_addtnl_wthnHH_trnsmssn_' num2str(j) '.csv']); % initial statuses of individuals from simulated data with ith presymptomatic/asymptomatic infectiousness
    hv=csvread(['PKDL_infctsnss' str '_addtnl_wthnHH_trnsmssn_' num2str(j) '.csv']); % infectiousness of PKDL cases from simulated data with ith presymptomatic/asymptomatic infectiousness
    nEmoves=round(sum(simdata.tI<108 & simdata.tR<=108)/5); % number of pre-symptomatic infection time updates per iteration
    nAupdts=round(sum(simdata.tA>0 & simdata.tA<=108)/5); % number of asymptomatic infection time updates per iteration
    VLStmMCMC2010FastLklhdMgrtnValidation(data,r1,p10,a,b,p2,uc{i},beta0,alpha0,epsilon0,delta0sc(i),lambda0,h0,h1,h2,h3,hmssng,h40,pI0,typ,niters,plotOutpt,rsltsc{i},para,simdata,Status0,hv,nEmoves,nAupdts)
end

%% CALCULATE DIC
delete(gcp('nocreate')) % delete current parallel worker pool
burnin=round(niters/5); % burn-in to discard
str1='DIC'; %str1='DIC_mean'; % string to prepend to rslts when saving DIC
ss='mode'; %ss='mean'; % summary statistic to use in calculation of DIC
RunCalcDIC(rsltsc,burnin,ss,str1)

%% CALCULATE DIC DIFFERENCES FROM BEST-FITTING MODEL
DICs=NaN(ndelta0s);
DICmin=NaN(1,ndelta0s);
DICdiffs=NaN(ndelta0s);
RMLs=NaN(ndelta0s);
DICminMdl=NaN(1,ndelta0s);
for i=1:ndelta0s
    [DICs(:,i),DICmin(i),DICdiffs(:,i),RMLs(:,i),DICminMdl(i)]=RunCalcDICdiffs(rsltsc(((i-1)*ndelta0s+1):i*ndelta0s),str1);
end

save('DICdiffs_addtnl_wthnHH_trnsmssn','DICdiffs')