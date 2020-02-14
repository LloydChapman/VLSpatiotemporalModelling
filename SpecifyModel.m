%SPECIFYMODEL Script to set input parameters for different models

%% DATA
% Load data
load('~/Dropbox/Visceral Leishmaniasis/CarynBernData/2010data/data_final2.mat'); % loads table called 'data'

% Set which paras are included in data
para=1:19;

%% MODEL PARAMETERS
% Parameters that can be estimated
% 1: beta - spatial transmission rate constant
% 2: alpha - distance scale factor for spatial kernel
% 3: epsilon - background transmission rate
% 4: delta - additional within-household transmission rate
% 5: lambda0 - historical asymptomatic infection rate
% 6: h4 - relative infectiousness of asymptomatics
% 7: pI - proportion of infections that lead to VL

% Set parameters for negative binomial incubation period distribution
r1=3; % shape parameter
mu=5; % guess for mean incubation period
p10=r1/(mu-1+r1); % initial guess for 'success' parameter
% Shape parameters for beta prior for p1
a=22;
b=(mu-1)*(a-1)/r1; % since mu-1=r1*b/(a-1)

% Set parameter for geometric asymptomatic infection period distribution
p2=1/5;

% Initial guesses/fixed values for transmission parameters
beta0=3;
alpha0=100;
epsilon0=1e-3;
delta0=1e-2;
pI0=0.15;
% Estimate lambda0 by fitting initial status model to LST data from 2002
s1=load('~/Dropbox/Visceral Leishmaniasis/CarynBernData/2004data/SpatiotemporalModelling/data_final.mat');
[pars,~]=FitCatModLST3(s1.data,p2);
lambda0=pars/12;

% Relative infectiousnesses of different forms of PKDL
h1=9/26/(10/15); % macular/papular
h3=18/21/(10/15); % nodular
h2=(h1+h3)/2; % plaque (assumed)
hmssng=(101/138*h1+31/138*h2+6/138*h3); % unexamined - use average relative infectiousness of PKDL cases for cases who weren't physically examined

%% MCMC parameters
niters=1e5; % number of iterations
plotOutpt=false; % flag for whether to plot output in real-time
runName='_AllParas';
 
%% PARAMETERS FOR DIFFERENT MODELS
% Names of models for saving
Mdls={'','_DELTA','_PRESX_NONINFCTS_ASX_NONINFCTS','_DELTA_PRESX_NONINFCTS_ASX_NONINFCTS','_HGHR_PRESX_ASX_INFCTS','_DELTA_HGHR_PRESX_ASX_INFCTS'};
if r1>1
    IPD='NBIP';
elseif r1==1
    IPD='GIP';
end
rslts=cellfun(@(x)['MCMC_' IPD '_PKDL_ASX' x runName],Mdls,'UniformOutput',false);
nMdls=numel(rslts);

% Set parameters to estimate in each model
u=repmat({1:3,1:4},1,nMdls/2);

% Set kernel type ('Cauchy' = Cauchy, 'Exp' = exponential, 'Const' = constant)
typ=repmat({'Exp'},1,nMdls);

% Set parameter values for the different models
beta0s=beta0*ones(1,nMdls);
alpha0s=alpha0*ones(1,nMdls);
delta0s=repmat([0,delta0],1,nMdls/2);
lambda0s=lambda0*ones(1,nMdls);
h0s=[0.01,0.01,0,0,0.02,0.02]; % relative infectiousness of pre-symptomatics
h1s=h1*ones(1,nMdls);
h2s=h2*ones(1,nMdls);
h3s=h3*ones(1,nMdls);
h40s=[0.01,0.01,0,0,0.02,0.02]; % relative infectiousness of asymptomatic individuals
pI0s=pI0*ones(1,nMdls);
hmssngs=hmssng*ones(1,nMdls);
inclLST=false(1,nMdls);
