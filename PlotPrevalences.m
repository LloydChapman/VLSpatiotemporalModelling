function PlotPrevalences(rslts,nsmpls,burnin1,varargin)

load(rslts)

if ~exist('z','var')
    z=burnin1+1:niters;
end

if nargin==3
    iters=randperm(numel(z),nsmpls);
else
    iters=varargin{1};
end

ziters=z(iters);

%% Construct matrices of event times that change across samples
% Asymptomatic infection
tAs1=tAs(:,ziters);

% Presymptomatic infection
tEs1=NaN(n,nsmpls);
tEs1(I,:)=tEs(:,ziters);

% VL onset
tIs=repmat(tI,[1,nsmpls]);
tIs(NONR,:)=tIsNONR(:,iters);
tIs(RNO,:)=tIsRNO(:,iters);
tIs(ANONR,:)=tIsANONR(:,iters);

% Recovery
% from asymptomatic infection
tRAs1=tRAs(:,ziters);

% from VL
tRs=repmat(tR,[1,nsmpls]);
tRs(NONR,:)=tRsNONR(:,ziters);
tRs(ONR,:)=tRsONR(:,ziters);
tRs(ANONR,:)=tRsANONR(:,ziters);
tRs(AONR,:)=tRsAONR(:,ziters);

% VL relapse and relapse treatment
tRLs=repmat(tRL,[1,nsmpls]);
tRLs(RLNO,:)=tRLsRLNO(:,ziters);
tRLRs=repmat(tRLR,[1,nsmpls]);
tRLRs(RLO,:)=tRLRsRLO(:,ziters);
tRLRs(RLNO,:)=tRLRsRLNO(:,ziters);

%% Count numbers in each state over time
numS=NaN(nsmpls,tmax);
numA=NaN(nsmpls,tmax);
numE=NaN(nsmpls,tmax);
numI=NaN(nsmpls,tmax);
numDor=NaN(nsmpls,tmax);
numP=NaN(nsmpls,tmax);
numR=NaN(nsmpls,tmax);
pop=NaN(nsmpls,tmax);
for j=1:tmax
    % Numbers in different states at time t
    numS(:,j)=sum(rng(:,1)-1<=j & rng(:,2)>j & tAs1>j & (tEs1>j|isnan(tEs1)) & ~prevK & ~actvK & ~IpreEXTIM & ~EXTIMsoonI & ~IpreINTIM & ~PpreINTIM & ~PpreEXTIM & ~(INTMIG_IN & isnan(tI) & isnan(tP) & tAs1==tmax+2),1);
    numA(:,j)=sum(tAs1<=j & rng(:,2)>j & tRAs1>j,1);
    numE(:,j)=sum(tEs1<=j & tIs>j,1);
    numI(:,j)=sum(((tIs<=j & tRs>j)|(tRLs<=j & (tRLRs>j|isnan(tRLRs)))) & rng(:,2)>j,1);
    numDor(:,j)=sum(tRs(IandP)<=j & tP(IandP)>j,1)+sum(tRs<=j & tRLs>j);
    numP(:,j)=sum(tP<=j & (tRP>j|isnan(tRP)));
    numR(:,j)=sum(tRs(setdiff(I,IandP))<=j & rng(setdiff(I,IandP),2)>j,1)+sum(tRP<=j & rng(:,2)>j,1)+sum(tRLRs<=j & rng(:,2)>j,1)+sum(tRAs1<=j & rng(:,2)>j,1);
    pop(:,j)=sum(rng(:,1)-1<=j & rng(:,2)>j);
end

%% Calculate prevalences
PrevS=numS./pop;
PrevA=numA./pop;
PrevE=numE./pop;
PrevI=numI./pop;
PrevDor=numDor./pop;
PrevP=numP./pop;
PrevR=numR./pop;

%% Plot numbers and prevalences of the different states over time
t=startyr+(0:tmax-1)/12;
% figure; hold on
% plot(t,numS')
% plot(t,numR')
% figure; hold on
% plot(t,numA')
% plot(t,numE')
% plot(t,numI')
% plot(t,numDor')
% plot(t,numP')
% figure; hold on 
% semilogy(t,numS')
% semilogy(t,numA')
% semilogy(t,numE')
% semilogy(t,numI')
% semilogy(t,numDor')
% semilogy(t,numP')
% semilogy(j,numR')
% figure; hold on
% plot(t,PrevS')
% plot(t,PrevR')
% figure; hold on
% plot(t,PrevA')
% plot(t,PrevE')
% plot(t,PrevI')
% plot(t,PrevDor')
% plot(t,PrevP')

% Define plot colours
clrs=[[254 224 77]/255;[254 195 87]/255;[245 150 79]/255;0.8 0.255 0.145;[173 163 198]/255;[81 130 187]/255;[146 208 88]/255];
nbins=50;

% Plot prevalence of S and R
x=NaN(nsmpls,tmax,2);
x(:,:,1)=PrevS;
x(:,:,2)=PrevR;
PlotCntrbtnModeAndHPDI(x,nbins,2,t,clrs([1,7],:),'Year','Prevalence',true,{'S','R'})
saveas(gcf,'PrevSandR')
saveas(gcf,'PrevSandR.eps','epsc')

% Plot prevalences of A, E, I, D, and P
x1=NaN(nsmpls,tmax,5);
x1(:,:,1)=PrevA;
x1(:,:,2)=PrevE;
x1(:,:,3)=PrevI;
x1(:,:,4)=PrevDor;
x1(:,:,5)=PrevP;
PlotCntrbtnModeAndHPDI(x1,nbins,5,t,clrs(2:6,:),'Year','Prevalence',true,{'A','E','I','D','P'})
saveas(gcf,'PrevAEIDP')
saveas(gcf,'PrevAEIDP.eps','epsc')
