clear

%% DATA
% Load data
load('C:\Users\timpo\OneDrive - University of Warwick\taubern_baybern\raw_data_plus_cleaning\matlab_bayesianmodel\data_final2.mat'); % loads table called 'data'
data.MOS_RX_NEW_SX = str2double(data.MOS_RX_NEW_SX);
% Set which paras are included in data
paras=1;

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
delta0=1e-2;
pI0=0.15;
npara=numel(paras);
% Estimate lambda0 by fitting initial status model to LST data from 2002
s1=load('C:\Users\timpo\OneDrive - University of Warwick\taubern_baybern\1999_2004dataset\data_final.mat');
[pars,~]=FitCatModLST3(s1.data,p2);
lambda0=pars/12;

h0=0.02; % relative infectiousness of pre-symptomatics
% Relative infectiousnesses of different forms of PKDL
h1=9/26/(10/15); % macular/papular
h3=18/21/(10/15); % nodular
h2=(h1+h3)/2; % plaque (assumed halfway between macular/papular and nodular)
hmssng=(101/138*h1+31/138*h2+6/138*h3); % unexamined - use average relative infectiousness of examined PKDL cases for cases who weren't physically examined
h40=0.02; % relative infectiousness of asymptomatics

% Set kernel type
typ='Exp';

%% MCMC PARAMETERS
niters=10; % number of iterations
plotOutpt=false; % flag for whether to plot output in real-time

tic
para = paras
runName=strcat('para',num2str(para));
IPD='NBIP';
rslts=['MCMC_' IPD '_PKDL_ASX' runName];
% Run MCMC
VLStmMCMC2010FastLklhdMgrtn(data,r1,p10,a,b,p2,u,beta0,alpha0,epsilon0,delta0,lambda0,h0,h1,h2,h3,hmssng,h40,pI0,typ,niters,plotOutpt,rslts,para)
toc
