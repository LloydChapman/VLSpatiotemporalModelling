function [M,HPDI,iters]=CalcCntrbtnMgrtnAsx(rslts,nsmpls,burnin1,varargin)

db='';
rng=[];

load(rslts)

% Load data
load('data_final2.mat')
% Select data for para
data=data(ismember(data.PARA,para),:);
% Rename longitude and latitude variables
data.Properties.VariableNames{'HHNEWLNG'}='longitude';
data.Properties.VariableNames{'HHNEWLAT'}='latitude';
% Calculate HH distance matrix
dHH=CalcHHDists(data);
dHHsqrd=[];

% Overwrite d0 in old output to make it work with new Knl_fast function
if size(d0,2)>size(d0,1)
    d0=speye(nHH);
end

if ~exist('z','var')
    z=burnin1+1:niters;
end

% Get indices of nsmpls to be drawn from the MCMC chain with the burn-in discarded
if nargin==3
    iters=randperm(numel(z),nsmpls);
else
    iters=varargin{1};
end

hP=zeros(nIPNIA,tmax);
for i=1:nIandP
    j=IandP(i);
    hP(I==j,max(tIM(j),tP(j))+1:min(min(tRP(j),tEM(j)),tmax))=hv(j);
end
for i=1:nPI
    j=PI(i);
    hP(nI+i,tP(j)+1:min(min(tRP(j),tEM(j)),tmax))=hv(j);
end
% Need to remove PKDL infectiousness for individual with simultaneous KA
% and PKDL here so we don't double count their contribution, as we treat 
% them as having KA infectiousness until they were treated for KA in MCMC
% code
PpreR=find(tR>tP&~isinf(tR));
hP(ismember(I,PpreR),tP(PpreR)+1:tR(PpreR))=0;
for i=1:nIMP
    j=IMP(i);
    hP(nI+nPI+nA+nIMI+i,tIM(j)+1:min(tRP(j),tmax))=hv(j);
end

% Create matrix for storing contribution from each infection state to total 
% force of infection (FOI) on all individuals
nsrcs=5; % number of infection sources (infectious states)
FOI=zeros(nsmpls,tmax,nsrcs);

% Create matrix for storing contribution from each infection state to total 
% FOI on all susceptible individuals
FOIonS=zeros(nsmpls,tmax,nsrcs);

% Create matrix for storing contribution from each infection state to FOIs 
% on KA cases at their infection times
nE=sum(sum(tEm)); % number of cases infected during study after initial incubation period window
FOIonE=zeros(nsmpls,nE,nsrcs);

I2=ismember(I,I1);
% figure(1);
for k=1:nsmpls
    m=z(iters(k));
    
    epsilon=p(m,3);
    
    tE=uint32(tEs(:,m)); % N.B. tEs has dimensions #KA cases x #iterations
    tEm=false(n,tmax);
    tEm((tE(I2)-1)*n+uint32(I1))=1;
    tA=uint32(tAs(:,m));
    tRA=uint32(tRAs(:,m));
    Asx=uint32(find(tA>=1 & tA<=tmax));
    tAm=false(n,tmax);
    tAm((tA(Asx)-1)*n+Asx)=1;
    prevA=uint32(find(tA==0 & tRA==0));
    actvA=uint32(find(tA==0 & tRA>0));
    IM_INprevAactvA=IM_IN(ismember(IM_OUT,[prevA;actvA]));
    
    S=1-max(preB,preIM)-max(max(max(cumsum(tEm,2),cumsum(tAm,2)),cumsum(tDm,2)),cumsum(tEMm,2));
    S(prevK,:)=0; % remove previous KA cases from susceptibles
    S(prevA,:)=0; % remove previously asymptomatically infected individuals from susceptibles
    S(actvA,:)=0; % remove initially actively asymptomatically infected individuals from susceptibles
    S(tA==0,:)=0; % remove previous and active asymptomatic infections from susceptibles
    S(tI<=maxIP,:)=0; % remove cases with onset before maxIP (also excludes cases with active KA at start of study and a couple of cases with onset in 2002 before immigration)
    S(IpreEXTIM,:)=0; % remove KA cases with onset before or at migration in
    S(EXTIMsoonI,:)=0; % remove KA cases with onset within 6 months of migration in
    S(IpreINTIM,:)=0; % remove KA cases with onset before internal migration in
    S(PpreINTIM,:)=0; % remove PKDL cases with onset before internal migration in 
    S(PpreEXTIM,:)=0; % remove PKDL cases (without prior KA) with onset before external migration in
    S(IM_INprevAactvA,:)=0; % remove susceptibility from 2nd observations of internal migrators asymptomatically infected before the start of the study
    S(IM_IN(ismember(IM_OUT,Asx)),:)=0; % remove susceptibility from 2nd observations of internal migrators asymptomatically infected during 1st observation

    tI(NONR)=tIsNONR(:,m);
    tRorD(NONR)=tRsNONR(:,m);
    tI(RNO)=tIsRNO(:,m);
    tRorD(ONR)=tRsONR(:,m);
    tI(ANONR)=tIsANONR(:,m);
    tRorD(ANONR)=tRsANONR(:,m);
    tRorD(AONR)=tRsAONR(:,m);
    tRLR(RLO)=tRLRsRLO(:,m);
    tRL(RLNO)=tRLsRLNO(:,m);
    tRLR(RLNO)=tRLRsRLNO(:,m);
    hE=zeros(nIPNIA,tmax);
    hI=zeros(nIPNIA,tmax);
    for i=1:nI
        j=I(i);
        hE(i,max(max(0,tIM(j)),tE(i))+1:tI(j))=h0;
        hI(i,max(tIM(j),tI(j))+1:min(tEM(j),tRorD(j)))=1;
    end
    for i=1:nA
        hI(nI+nPI+i,1:tRorD(A(i)))=1;
    end
    for i=1:nIMI
        j=IMI(i);
        hI(nI+nPI+nA+i,tIM(j)+1:tRorD(j))=1;
    end
    for i=1:nRL
        j=RL(i);
        hI(IPNIA==j,tRL(j)+1:min(min(tRLR(j),tEM(j)),tmax))=1;
    end
    hA=zeros(n,tmax);
    A1=find(tA>=0 & tA<tmax+1 & tRA>0 & isnan(tP));
    for i=1:numel(A1)
        j=A1(i);
        hA(j,tA(j)+1:min(min(min(tRA(j),tEM(j)),tD(j)),tmax))=p(m,6);
    end
    RAobs1actvA=find(tA==0 & tRA>=tEM & INTMIG_OUT);
    RAobs2actvA=IM_IN(ismember(IM_OUT,RAobs1actvA));    
    RAobs1=find(tA>0 & tA<tEM & tRA>=tEM & INTMIG_OUT);
    RAobs2=IM_IN(ismember(IM_OUT,RAobs1));
    for i=1:numel(RAobs2actvA)
        j=RAobs1actvA(i);
        j1=RAobs2actvA(i);
        hA(j1,rng(j,2)+1:min(min(tRA(j),rng(j1,2)),tmax))=p(m,6);
    end
    for i=1:numel(RAobs2)
        j=RAobs1(i);
        j1=RAobs2(i);
        hA(j1,rng(j,2)+1:min(min(tRA(j),rng(j1,2)),tmax))=p(m,6);
    end
    hA=sparse(hA);
    hPA=zeros(nPA,tmax);
    for i=1:nPA
        j=PA(i);
        hPA(i,max(1,tA(j)+1):tRA(j))=p(m,6);
        hPA(i,tP(j)+1:min(min(tRP(j),tEM(j)),tmax))=hv(j);
    end
    
    KHH=Knl_fast(dHH,dHHsqrd,p(m,2),p(m,1),typ,n,nHH,ib,f,1:n);
    rateHHA=p(m,1)*KHH;
    if ismember(4,u)
        rateHHA=rateHHA+p(m,4)*d0;
    end
    rateHH=rateHHA(:,ib(IPNIA));
    rateHHPA=rateHHA(:,ib(PA));
    lambdaHHA=rateHHA(:,ib)*hA;
    lambdaHHE=rateHH*hE;
    lambdaHHI=rateHH*hI;
    lambdaHHP=rateHH*hP+rateHHPA*hPA;
    lambdaA=lambdaHHA(ib,:);
    lambdaE=lambdaHHE(ib,:);
    lambdaI=lambdaHHI(ib,:);
    lambdaP=lambdaHHP(ib,:);
    lambdaAE=lambdaA(tEm);
    lambdaEE=lambdaE(tEm);
    lambdaIE=lambdaI(tEm);
    lambdaPE=lambdaP(tEm);
    
    % Calculate contributions from background, asymptomatics, 
    % presymptomatics, KA cases and PKDL cases to:
    
    % total FOI on all individuals
    FOI(k,:,1)=n*epsilon;
    FOI(k,:,2)=sum(lambdaA);
    FOI(k,:,3)=sum(lambdaE);
    FOI(k,:,4)=sum(lambdaI);    
    FOI(k,:,5)=sum(lambdaP);
    
    % total FOI on all susceptible individuals 
    FOIonS(k,:,1)=sum(epsilon*S);
    FOIonS(k,:,2)=sum(lambdaA.*S);
    FOIonS(k,:,3)=sum(lambdaE.*S);
    FOIonS(k,:,4)=sum(lambdaI.*S);
    FOIonS(k,:,5)=sum(lambdaP.*S);
    
    % FOIs on KA cases at their infection times 
    FOIonE(k,:,1)=epsilon;
    FOIonE(k,:,2)=lambdaAE;
    FOIonE(k,:,3)=lambdaEE;
    FOIonE(k,:,4)=lambdaIE;
    FOIonE(k,:,5)=lambdaPE;
    
end

%% PLOTS
% Define plot colours
clrs=[0.96 0.9 0.8;[254 195 87]/255;[245 150 79]/255;0.8 0.255 0.145;[81 130 187]/255];
nclrs=size(clrs,1);
% Make time vector (in years)
t=startyr+(0:tmax-1)/12;

% figure; plot(t,mean(Sus))
% xlabel('Time'); ylabel('No. susceptibles')
% saveas(gcf,'NumSscptblesOverTime.png','png')

nbins=50;
lgd={'Background','Asx','Presx','VL','PKDL'};

%% Plot FOI on whole population
PlotCntrbtnModeAndHPDI(FOI,nbins,nsrcs,t,clrs,'Year','Total infection pressure (month^{-1})',true,lgd)
saveas(gcf,'AbsltCntrbtnFOIAsx')
saveas(gcf,'AbsltCntrbtnFOIAsx.eps','epsc')
saveas(gcf,'AbsltCntrbtnFOIAsx.png')
PlotCntrbtnModeAndHPDI(bsxfun(@rdivide,FOI,sum(FOI,3)),nbins,nsrcs,t,clrs,'Year','Relative contribution to total infection pressure',true,lgd,'west')
saveas(gcf,'RltveCntrbtnFOIAsx')
saveas(gcf,'RltveCntrbtnFOIAsx.eps','epsc')
saveas(gcf,'RltveCntrbtnFOIAsx.png')

%% Plot FOI on susceptibles
PlotCntrbtnModeAndHPDI(FOIonS,nbins,nsrcs,t,clrs,'Year','Total infection pressure (month^{-1})',true,lgd)
saveas(gcf,'AbsltCntrbtnFOIonSAsx')
saveas(gcf,'AbsltCntrbtnFOIonSAsx.eps','epsc')
saveas(gcf,'AbsltCntrbtnFOIonSAsx.png')
PlotCntrbtnModeAndHPDI(bsxfun(@rdivide,FOIonS,sum(FOIonS,3)),nbins,nsrcs,t,clrs,'Year','Relative contribution to total infection pressure',true,lgd,'west')
saveas(gcf,'RltveCntrbtnFOIonSAsx')
saveas(gcf,'RltveCntrbtnFOIonSAsx.eps','epsc')
saveas(gcf,'RltveCntrbtnFOIonSAsx.png')

%% Plot FOIs on KA cases at infection times
totFOIonE=sum(FOIonE,3);
PlotCntrbtnModeAndHPDI(FOIonE,nbins,nsrcs,1:nE,clrs,'Case number (by onset)','Infection pressure on case (month^{-1})',false,lgd)
saveas(gcf,'AbsltIndvdlCntrbtnAsx')
saveas(gcf,'AbsltIndvdlCntrbtnAsx.eps','epsc')
saveas(gcf,'AbsltIndvdlCntrbtnAsx.png')
PlotCntrbtnModeAndHPDI(bsxfun(@rdivide,FOIonE,totFOIonE),nbins,nsrcs,1:nE,clrs,'Case number (by onset)','Relative contribution to infection pressure',false,lgd,'west')
saveas(gcf,'RltveIndvdlCntrbtnAsx')
saveas(gcf,'RltveIndvdlCntrbtnAsx.eps','epsc')
saveas(gcf,'RltveIndvdlCntrbtnAsx.png')

%% Overall contribution to new KA cases over the course of the epidemic
figure; hold on
M=zeros(nsrcs,1);
HPDI=zeros(nsrcs,2);
for i=1:nsrcs
    hf6=histogram(sum(FOIonE(:,:,i),2)./sum(totFOIonE,2),nbins,'Normalization','pdf');
    hf6.LineStyle='none';
    hf6.FaceColor=clrs(i,:);
    hf6.FaceAlpha=0.8;
    [M(i),HPDI(i,:)]=CalcModeAndHPDI(sum(FOIonE(:,:,i),2)./sum(totFOIonE,2),nbins);
end
set(gca,'FontSize',14)
xlabel('Relative overall contribution')
ylabel('Density')
legend('Bckgrnd','Asx','Presx','VL','PKDL')
hold off
saveas(gcf,'RltveOvrlCntrbtnAsx')
saveas(gcf,'RltveOvrlCntrbtnAsx.eps','epsc')
saveas(gcf,'RltveOvrlCntrbtnAsx.png')

%% Save mode and HPDI for relative overall contribution of each infection state to FOIs on KA cases at their infection times
save('ModeAndHPDIRltveOvrlCntrbtnAsx','M','HPDI','iters')
