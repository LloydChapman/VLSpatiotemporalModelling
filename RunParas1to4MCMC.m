%RUNPARA1MCMC Script to run MCMC on para 1 data with different asymptomatic
%proportions for use in validation
clear

%% DATA
% Load data
load('data_final2.mat'); % loads table called 'data'

% Set which paras are included in data
para=1:4;

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
pI0s=[0.1,0.15,0.2];
npI0s=numel(pI0s);
% Estimate lambda0 by fitting initial status model to LST data from 2002
s1=load('data_final');
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
niters=1e5; % number of iterations
plotOutpt=false; % flag for whether to plot output in real-time
runName='_Paras1to4';
% Names of runs for saving
if r1>1
    IPD='NBIP';
elseif r1==1
    IPD='GIP';
end
rslts=['MCMC_' IPD '_PKDL_ASX' runName];

%% RUN MCMC ON OBSERVED DATA WITH DIFFERENT PROPORTIONS OF ASYMPTOMATIC INFECTION
% Run MCMC chains in parallel
parfor i=1:npI0s
% for i=1:npI0s
    % Set maximum number of computational threads for each parallel worker
    maxNumCompThreads(floor(12/npI0s))
    % Run MCMC
    VLStmMCMC2010FastLklhdMgrtn(data,r1,p10,a,b,p2,u,beta0,alpha0,epsilon0,delta0,lambda0,h0,h1,h2,h3,hmssng,h40,pI0s(i),typ,niters,plotOutpt,[rslts '_' num2str(i)],para)
end

%% SAVE PARAMETER VALUES AND MISSING DATA FOR SAMPLE CLOSE TO POSTERIOR MEAN
for i=1:npI0s
    % Load results
    load([rslts '_' num2str(i)],'p','p1','burnin','niters')
    % Vector for removing burn-in
    z=burnin+1:niters;
    % Calculate marginal posterior means
    mean_p=mean(p(z,:),1);
    mean_p1=mean(p1(z));
    % Find sample closest to posterior means
    x=sum(abs(bsxfun(@minus,[p(z,:),p1(z)],[mean_p,mean_p1])./[mean_p,mean_p1]),2);
    iter=find(x==min(x));
    % Save parameter values and missing data for this iteration
    SaveSamplesAsCSVs([rslts '_' num2str(i)],[],burnin,[runName '_' num2str(i)],iter)
end