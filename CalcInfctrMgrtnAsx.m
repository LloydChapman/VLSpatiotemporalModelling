function [infctn,infctr,src,dists,times,SI,onsetinfctr,rcvryinfctr,meandists,meantimes,infctrmax,srcmax,distsmax,timesmax,SImax,onsetinfctrmax,rcvryinfctrmax,meandistsmax,meantimesmax,RjA_I,RjI_I,RjP_I,Rj_I,Rts_I,Rt_I,diA_I,diI_I,diP_I,di_I,tiA_I,tiI_I,tiP_I,ti_I,infctnA,infctrA,srcA,distsA,timesA,SIA,onsetinfctrA,rcvryinfctrA,infctrmaxA,srcmaxA,distsmaxA,timesmaxA,SImaxA,onsetinfctrmaxA,rcvryinfctrmaxA,RjA_A,RjI_A,RjP_A,Rj_A,Rts_A,Rt_A,diA_A,diI_A,diP_A,di_A,tiA_A,tiI_A,tiP_A,ti_A,iters,OT,Rj,Rts,Rt]=CalcInfctrMgrtnAsx(rslts,nsmpls,burnin1,varargin)

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

if nargin==3
    iters=randperm(numel(z),nsmpls);
else
    iters=varargin{1};
end

% Construct infectiousness and status matrices
% States in status matrix numbered as follows:
% 0=unborn/dead
% 1=susceptible
% 2=asymptomatic
% 3=pre-symptomatic
% 4=KA
% 5=dormant
% 6=PKDL
% 7=recovered
% 9=background
stat=[zeros(n,1),rngm(:,1:end-1)];
stat(Pres0,1)=1;
h=zeros(nIPNIA,tmax);
for i=1:nIMI
    j=IMI(i);
    h(nI+nPI+nA+i,tIM(j)+1:tRorD(j))=1;
    stat(j,tIM(j)+1:tRorD(j))=4;
    stat(j,tRorD(j)+1:tmax)=7;
end
% Fill in status matrix for individuals who had PKDL without prior KA or
% with onset before startyr as these will not change with
% infection/KA onset/recovery time updates
hP=zeros(nIPNIA,tmax);
for i=1:nIandP
    j=IandP(i);
    hP(I==j,max(tIM(j),tP(j))+1:min(min(tRP(j),tEM(j)),tmax))=hv(j);
    stat(j,max(tIM(j),tP(j))+1:min(min(tRP(j),tEM(j)),tmax))=6;
    stat(j,min(tRP(j)+1,tmax+1):tmax)=7;
end
for i=1:nPI
    j=PI(i);
    hP(nI+i,tP(j)+1:min(min(tRP(j),tEM(j)),tmax))=hv(j);
    if ismember(j,find(prevK|actvK))
        stat(j,1:tP(j))=5;
    end
    stat(j,tP(j)+1:min(min(tRP(j),tEM(j)),tmax))=6;
    stat(j,min(tRP(j)+1,tmax+1):tmax)=7;
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
    stat(j,tIM(j)+1:tRP(j))=6;
    stat(j,tRP(j)+1:tmax)=7;
end

% Make index vector of cases with onset after initial window, i.e. whose
% infection sources are imputed
nI1=numel(I1);
infctn=NaN(nI1,nsmpls);
infctr=NaN(nI1,nsmpls);
src=NaN(nI1,nsmpls);
infctrmax=NaN(nI1,nsmpls);
srcmax=NaN(nI1,nsmpls);
RjA_I=NaN(nsmpls,n);
RjI_I=NaN(nsmpls,nIPNIA);
RjP_I=NaN(nsmpls,nIPNIA+nPA);
Rj_I=NaN(nsmpls,n+nIPNIA+nPA+1);
Rts_I=NaN(3,tmax-maxIP,nsmpls);
Rt_I=zeros(nsmpls,tmax-maxIP);
diA_I=NaN(nI1,nsmpls);
diI_I=NaN(nI1,nsmpls);
diP_I=NaN(nI1,nsmpls);
di_I=NaN(nI1,nsmpls);
tiA_I=NaN(nI1,nsmpls);
tiI_I=NaN(nI1,nsmpls);
tiP_I=NaN(nI1,nsmpls);
ti_I=NaN(nI1,nsmpls);
dists=NaN(nI1,nsmpls);
onsetinfctr=NaN(nI1,nsmpls);
rcvryinfctr=NaN(nI1,nsmpls);
times=NaN(nI1,nsmpls);
SI=NaN(nI1,nsmpls);
distsmax=NaN(nI1,nsmpls);
onsetinfctrmax=NaN(nI1,nsmpls);
rcvryinfctrmax=NaN(nI1,nsmpls);
timesmax=NaN(nI1,nsmpls);
SImax=NaN(nI1,nsmpls);
I2=ismember(I,I1);
infctnA=NaN(n,nsmpls);
infctrA=NaN(n,nsmpls);
srcA=NaN(n,nsmpls);
infctrmaxA=NaN(n,nsmpls);
srcmaxA=NaN(n,nsmpls);
RjA_A=NaN(nsmpls,n);
RjI_A=NaN(nsmpls,nIPNIA);
RjP_A=NaN(nsmpls,nIPNIA+nPA);
Rj_A=NaN(nsmpls,n+nIPNIA+nPA+1);
Rts_A=NaN(3,tmax-maxIP,nsmpls);
Rt_A=zeros(nsmpls,tmax-maxIP);
diA_A=NaN(n,nsmpls);
diI_A=NaN(n,nsmpls);
diP_A=NaN(n,nsmpls);
di_A=NaN(n,nsmpls);
tiA_A=NaN(n,nsmpls);
tiI_A=NaN(n,nsmpls);
tiP_A=NaN(n,nsmpls);
ti_A=NaN(n,nsmpls);
distsA=NaN(n,nsmpls);
onsetinfctrA=NaN(n,nsmpls);
rcvryinfctrA=NaN(n,nsmpls);
timesA=NaN(n,nsmpls);
SIA=NaN(n,nsmpls);
distsmaxA=NaN(n,nsmpls);
onsetinfctrmaxA=NaN(n,nsmpls);
rcvryinfctrmaxA=NaN(n,nsmpls);
timesmaxA=NaN(n,nsmpls);
SImaxA=NaN(n,nsmpls);
OT=NaN(nsmpls,n+2*nIPNIA+nPA+1);
Rt=zeros(nsmpls,tmax-maxIP);
Rts=zeros(3,tmax-maxIP,nsmpls);
for k=1:nsmpls
    m=z(iters(k));
    % Get current infection, onset and recovery times
    tE=double(tEs(:,m));
    tA=double(tAs(:,m));
    tA(tA==tmax+2)=NaN;
    tRA=double(tRAs(:,m));
    tRA(tRA==tmax+2)=NaN;
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
    
    A1=find(tA>=0 & tA<tmax+1 & tRA>0 & isnan(tP));
    RAobs1actvA=find(tA==0 & tRA>=tEM & INTMIG_OUT);
    RAobs2actvA=IM_IN(ismember(IM_OUT,RAobs1actvA));    
    RAobs1=find(tA>0 & tA<tEM & tRA>=tEM & INTMIG_OUT);
    RAobs2=IM_IN(ismember(IM_OUT,RAobs1));
    % Overwrite rows corresponding to KA cases (with onset during study or 
    % active KA at start of study) 
    h([1:nI,nI+nPI+1:nI+nPI+nA],:)=zeros(nI+nA,tmax);
    stat([I;A;A1;RAobs2actvA;RAobs2;PA],:)=[zeros(numel([I;A;A1;RAobs2actvA;RAobs2;PA]),1),rngm([I;A;A1;RAobs2actvA;RAobs2;PA],1:end-1)];
    stat(intersect([I;A;A1;RAobs2actvA;RAobs2;PA],Pres0),1)=1;
    % Fill in current infectiousness and status matrices
    for i=1:nI
        j=I(i);
        h(i,max(max(0,tIM(j)),tE(i))+1:tI(j))=h0;
        h(i,max(tIM(j),tI(j))+1:min(tEM(j),tRorD(j)))=1;
        stat(j,max(max(0,tIM(j)),tE(i))+1:tI(j))=3;
        stat(j,max(tIM(j),tI(j))+1:min(tEM(j),tRorD(j)))=4;
        if ismember(j,IandP)
            stat(j,max(tIM(j),tRorD(j))+1:tP(j))=5;
            stat(j,max(tIM(j),tP(j))+1:min(min(tRP(j),tEM(j)),tmax))=6;
            stat(j,min(tRP(j)+1,tmax+1):tmax)=7;
        else
            stat(j,tRorD(j)+1:min(min(tEM(j),tD(j)),tmax))=7;
        end
    end
    for i=1:nA
        j=A(i);
        h(nI+nPI+i,1:tRorD(j))=1;
        stat(j,1:tRorD(j))=4;
        % status of cases with later PKDL gets overwritten above, so need to reinstate it here
        if ismember(j,intersect(PI,A)) % case had later PKDL
           stat(j,max(1,tRorD(j)+1):tP(j))=5;
           stat(j,tP(j)+1:min(min(tRP(j),tEM(j)),tmax))=6;
           stat(j,min(tRP(j)+1,tmax+1):tmax)=7;
        else % case didn't have later PKDL
           stat(j,max(1,tRorD(j)+1):min(min(tEM(j),tD(j)),tmax))=7;
        end
    end
    for i=1:nRL
       j=RL(i);
       h(I==j,tRL(j)+1:min(min(tRLR(j),tEM(j)),tmax))=1;
       stat(j,tRL(j)+1:min(min(tRLR(j),tEM(j)),tmax))=4;
    end
    % Completely overwrite asymptomatic infectiousness matrix
    hA=zeros(n,tmax);    
    for i=1:numel(A1)
        j=A1(i);
        hA(j,tA(j)+1:min(min(min(tRA(j),tEM(j)),tD(j)),tmax))=p(m,6);
        stat(j,tA(j)+1:min(min(min(tRA(j),tEM(j)),tD(j)),tmax))=2;
        stat(j,tRA(j)+1:min(min(tEM(j),tD(j)),tmax))=7;
    end
    for i=1:numel(RAobs2actvA)
        j=RAobs1actvA(i);
        j1=RAobs2actvA(i);
        hA(j1,rng(j,2)+1:min(min(tRA(j),rng(j1,2)),tmax))=p(m,6);
        stat(j1,rng(j,2)+1:min(min(tRA(j),rng(j1,2)),tmax))=2;
        stat(j1,tRA(j)+1:min(rng(j1,2),tmax))=7;
    end
    for i=1:numel(RAobs2)
        j=RAobs1(i);
        j1=RAobs2(i);
        hA(j1,rng(j,2)+1:min(min(tRA(j),rng(j1,2)),tmax))=p(m,6);
        stat(j1,rng(j,2)+1:min(min(tRA(j),rng(j1,2)),tmax))=2;
        stat(j1,tRA(j)+1:min(rng(j1,2),tmax))=7;
    end
    hA=sparse(hA);
    hPA=zeros(nPA,tmax);
    for i=1:nPA
        j=PA(i);
        hPA(i,max(1,tA(j)+1):tRA(j))=p(m,6);
        hPA(i,tP(j)+1:min(min(tRP(j),tEM(j)),tmax))=hv(j);
        stat(j,max(1,tA(j)+1):tRA(j))=2;
        stat(j,max(1,tRA(j)+1):tP(j))=5;
        stat(j,tP(j)+1:min(min(tRP(j),tEM(j)),tmax))=6;
        stat(j,min(min(tRP(j),tEM(j)),tmax)+1:min(min(tEM(j),tD(j)),tmax))=7;
    end
    
    % Store pre-symptomatic infection times
    infctn(:,k)=tE(I2);
    Asx=find(tA>=1 & tA<=tmax);
    nAsx=numel(Asx);
    infctnA(:,k)=tA;
    
    % Calculate pairwise infectious pressures at infection times
    KHH=Knl_fast(dHH,dHHsqrd,p(m,2),p(m,1),typ,n,nHH,ib,f,1:n);
    rateHHA=p(m,1)*KHH;
    if ismember(4,u)
        rateHHA=rateHHA+p(m,4)*d0;
    end
    tmp=rateHHA(ib,:);
    rate=tmp(:,ib(IPNIA));
    rateA=tmp(:,ib);
    ratePA=tmp(:,ib(PA));
    
    % Make vector of "onset" times of infectious individuals to pass into function to calculate reproduction numbers
    onsetI=[tI(I);tP(PI);tI(A);tIM(IMI);tIM(IMP)]; % use immigration time for "onset" time of individuals that had KA/PKDL at time of immigration
    onsetP=[tP(I);tP(PI);tP(A);tIM(IMI);tIM(IMP)]; % use immigration time for "onset" time of individuals that had KA/PKDL at time of immigration   
    
    [infctr(:,k),src(:,k),infctrmax(:,k),srcmax(:,k),RjA_I(k,:),RjI_I(k,:),RjP_I(k,:),Rj_I(k,:),Rts_I(:,:,k),Rt_I(k,:),diA_I(:,k),diI_I(:,k),diP_I(:,k),di_I(:,k),tiA_I(:,k),tiI_I(:,k),tiP_I(:,k),ti_I(:,k)]=GetInfctrAndSrc(infctn(:,k),rateA,rate,ratePA,hA,h,hP,hPA,p(m,:),I1,[A1;RAobs2actvA;RAobs2],nI1,tA,onsetI,onsetP,tP(PA),dHH,ib,PA,IPNIA,maxIP,tmax,stat,n,nIPNIA);
    [infctrA(Asx,k),srcA(Asx,k),infctrmaxA(Asx,k),srcmaxA(Asx,k),RjA_A(k,:),RjI_A(k,:),RjP_A(k,:),Rj_A(k,:),Rts_A(:,:,k),Rt_A(k,:),diA_A(Asx,k),diI_A(Asx,k),diP_A(Asx,k),di_A(Asx,k),tiA_A(Asx,k),tiI_A(Asx,k),tiP_A(Asx,k),ti_A(Asx,k)]=GetInfctrAndSrc(infctnA(Asx,k),rateA,rate,ratePA,hA,h,hP,hPA,p(m,:),Asx,[A1;RAobs2actvA;RAobs2],nAsx,tA,onsetI,onsetP,tP(PA),dHH,ib,PA,IPNIA,maxIP,tmax,stat,n,nIPNIA);

    % Calculate time-dependent population-level effective reproduction
    % numbers - this part must be inside loop as onset changes in each
    % iteration
    onset=[tA;onsetI;tP(PA);NaN];
    [Rt(k,:),Rts(:,:,k)]=CalcRt(Rj_I(k,:)+Rj_A(k,:),RjA_I(k,:)+RjA_A(k,:),RjI_I(k,:)+RjI_A(k,:),RjP_I(k,:)+RjP_A(k,:),onset,tA,onsetI,maxIP,tmax);
    
    OT(k,:)=[tRA-tA;min(tRorD(I),tEM(I))-tI(I);NaN(nPI,1);min(tRorD(A),tEM(A))-tI(A);min(tRorD(IMI),tEM(IMI))-tIM(IMI);NaN(nIMP,1);min(min(tRP(I),tEM(I)),tmax)-tP(I);min(min(tRP(PI),tEM(PI)),tmax)-tP(PI);NaN(nA,1);NaN(nIMI,1);min(min(tRP(IMP),tEM(IMP)),tmax)-tIM(IMP);min(min(tRP(PA),tEM(PA)),tmax)-tP(PA);NaN]';
    % deal with internal migrators who had PKDL in their 2nd observation 
    % having PKDL onset and treatment times in their 1st observations - a 
    % better fix would be to restrict the part for onset times of PKDL 
    % cases to just PKDL cases, and likewise for the part for VL onset times
    OT(k,n+nIPNIA+find(~ismember(I,IandP)))=NaN;
    
    % Get infection distances and intervals for infectees in each tree
    [dists(:,k),onsetinfctr(:,k),rcvryinfctr(:,k),times(:,k),SI(:,k)]=GetDistsAndTimes(infctn(:,k),infctr(:,k),src(:,k),dHH,ib,I1,tA,tI,tP,tRorD,tRP);
    [distsA(Asx,k),onsetinfctrA(Asx,k),rcvryinfctrA(Asx,k),timesA(Asx,k),SIA(Asx,k)]=GetDistsAndTimes(infctnA(Asx,k),infctrA(Asx,k),srcA(Asx,k),dHH,ib,Asx,tA,tI,tP,tRorD,tRP);

    % Get infection distances and intervals for most likely infectors in each tree
    [distsmax(:,k),onsetinfctrmax(:,k),rcvryinfctrmax(:,k),timesmax(:,k),SImax(:,k)]=GetDistsAndTimes(infctn(:,k),infctrmax(:,k),srcmax(:,k),dHH,ib,I1,tA,tI,tP,tRorD,tRP);
    [distsmaxA(Asx,k),onsetinfctrmaxA(Asx,k),rcvryinfctrmaxA(Asx,k),timesmaxA(Asx,k),SImaxA(Asx,k)]=GetDistsAndTimes(infctnA(Asx,k),infctrmaxA(Asx,k),srcmaxA(Asx,k),dHH,ib,Asx,tA,tI,tP,tRorD,tRP);    
    
end

%% Calculate total individual-level effective reproduction numbers (sum of R's for new VL cases and new asx infctns)
Rj=Rj_I+Rj_A;
idxI=[find(ismember(I,I1));nI+nPI+nA+1:nI+nPI+nA+nIMI]; % include IMI as their onset time > maxIP
RjI_I=RjI_I(:,idxI);
RjI_A=RjI_A(:,idxI);
RjI=RjI_I+RjI_A;
idxP=find(ismember([I;PI;NaN(nA,1);IMI;IMP;PA],[IandP;PI;IMP;PA])); % replace indices for initially active KA cases by NaN to avoid duplication of initially active KA cases that later developed PKDL
RjP_I=RjP_I(:,idxP);
RjP_A=RjP_A(:,idxP);
RjP=RjP_I+RjP_A;
RjA=RjA_I+RjA_A;

%% Calculate distances to infectors for VL cases & asx infctns
di=[di_I;di_A];
diA=[diA_I;diA_A];
diI=[diI_I;diI_A];
diP=[diP_I;diP_A];
ti=[ti_I;ti_A];
tiA=[tiA_I;tiA_A];
tiI=[tiI_I;tiI_A];
tiP=[tiP_I;tiP_A];

%% Set colours for plotting
dfltclrs=get(gca,'ColorOrder');
clrs=[0.96 0.9 0.8;[254 195 87]/255;[245 150 79]/255;0.8 0.255 0.145;[81 130 187]/255];
clrs=clrs([2,4,5],:);

%% Plot individual-level effective reproduction numbers for VL and PKDL cases
figure;
edges=0:ceil(max(max(mean(RjI,1)),max(mean(RjP,1))));
PlotRjDistn(RjI,1,'Mean no. secondary infections','tex',clrs(2,:),edges); hold on
PlotRjDistn(RjP,1,'Mean no. secondary infections','tex',clrs(3,:),edges) 
legend(gca,{'VL','PKDL'})
hold off
saveas(gcf,'RjDistn')
saveas(gcf,'RjDistn.eps','epsc')
saveas(gcf,'RjDistn.png')
meanRjI=mean(RjI,1);
minRjI=min(meanRjI);
maxRjI=max(meanRjI);
[parsI,CII]=FitGammaDistnToRjDistn(meanRjI);
meanRjP=mean(RjP,1);
minRjP=min(meanRjP);
maxRjP=max(meanRjP);
[parsP,CIP]=FitGammaDistnToRjDistn(meanRjP);

[~,iI]=sort(tI([I1;IMI]));
PlotRj(RjI(:,iI),1:numel(I1)+nIMI,'VL case number (by onset)','No. secondary infections')
set(gcf,'Position',[0 50 round(scrnsz(3)) round(scrnsz(3)/2)]);
saveas(gcf,'RjI')
saveas(gcf,'RjI.eps','epsc')
saveas(gcf,'RjI.png')

[~,iP]=sort(tP([IandP;PI;IMP;PA]));
PlotRj(RjP(:,iP),1:nIandP+nPI+nIMP+nPA,'PKDL case number (by onset)','No. secondary infections')
set(gcf,'Position',[0 50 round(scrnsz(3)) round(scrnsz(3)/2)]);
saveas(gcf,'RjP')
saveas(gcf,'RjP.eps','epsc')
saveas(gcf,'RjP.png')

%% Plot individual-level effective reproduction numbers for asymptomatic individuals
figure;
PlotRjDistn(RjA,1,'Mean no. secondary infections','tex',clrs(1,:))
saveas(gcf,'RjADistn')
saveas(gcf,'RjADistn.eps','epsc')
saveas(gcf,'RjADistn.png')
meanRjA=mean(RjA,1,'omitnan');
min(meanRjA);
max(meanRjA);

%% Plot effective reproduction numbers against onset-to-recovery time
figure;
plot(mean(OT(:,1:n),1,'omitnan'),meanRjA,'.','Color',clrs(1,:),'MarkerSize',14); hold on
plot(mean(OT(:,n+idxI),1),meanRjI,'.','Color',clrs(2,:),'MarkerSize',14); hold on
plot(mean(OT(:,n+nIPNIA+idxP),1),meanRjP,'.','Color',clrs(3,:),'MarkerSize',14)
set(gca,'FontSize',14)
legend(gca,{'Asx','VL','PKDL'},'Location','NorthWest')
xlabel('Onset-to-recovery time (months)'); ylabel('Mean no. secondary infections')
hold off
saveas(gcf,'RjVsOT')
saveas(gcf,'RjVsOT.eps','epsc')
saveas(gcf,'RjVsOT.png')

%% Plot VL and PKDL case reproduction numbers for new infections leading to VL
figure;
edges1=0:0.2:ceil(max(max(mean(RjI_I,1)),max(mean(RjP_I,1))));
PlotRjDistn(RjI_I,1,'Mean no. secondary VL cases','tex',clrs(2,:),edges1); hold on
PlotRjDistn(RjP_I,1,'Mean no. secondary VL cases','tex',clrs(3,:),edges1);
legend(gca,{'VL','PKDL'})
hold off
saveas(gcf,'RjVLDistn')
saveas(gcf,'RjVLDistn.eps','epsc')
saveas(gcf,'RjVLDistn.png')
meanRjI_I=mean(RjI_I,1);
minRjI_I=min(meanRjI_I);
maxRjI_I=max(meanRjI_I);
[parsI_I,CII_I]=FitGammaDistnToRjDistn(meanRjI_I);
meanRjP_I=mean(RjP_I,1);
minRjP_I=min(meanRjP_I);
maxRjP_I=max(meanRjP_I);
[parsP_I,CIP_I]=FitGammaDistnToRjDistn(meanRjP_I);
PlotRj(RjI_I(:,iI),1:numel(I1)+nIMI,'VL case number (by onset)','No. secondary VL cases')
saveas(gcf,'RjIVL')
saveas(gcf,'RjIVL.eps','epsc')
saveas(gcf,'RjIVL.png')
PlotRj(RjP_I(:,iP),1:nIandP+nPI+nIMP+nPA,'PKDL case number (by onset)','No. secondary VL cases')
saveas(gcf,'RjPVL')
saveas(gcf,'RjPVL.eps','epsc')
saveas(gcf,'RjPVL.png')

%% Plot median effective reproduction number over time
figure;
t=startyr+(maxIP:tmax-1)/12;
hf2(4)=PlotMCTAndCI(t,median(Rt,1),[quantile(Rt,0.025,1);quantile(Rt,0.975,1)],[0.5 0.5 0.5],0.2,'Time','');
for i=1:3
    hold on;
    hf2(i)=PlotMCTAndCI(t,median(Rts(i,:,:),3),[quantile(Rts(i,:,:),0.025,3);quantile(Rts(i,:,:),0.975,3)],clrs(i,:),0.2,'Time','$$R(t)$$');
end
legend(hf2,{'Asx','VL','PKDL','All'})
saveas(gcf,'Rt')
saveas(gcf,'Rt.eps','epsc')
saveas(gcf,'Rt.png')

%% Plot median "effective reproduction number" for new infections leading to VL
figure;
hf3(4)=PlotMCTAndCI(t,median(Rt_I,1),[quantile(Rt_I,0.025,1);quantile(Rt_I,0.975,1)],[0.5 0.5 0.5],0.2,'Time','');
for i=1:3
    hold on;
    hf3(i)=PlotMCTAndCI(t,median(Rts_I(i,:,:),3),[quantile(Rts_I(i,:,:),0.025,3);quantile(Rts_I(i,:,:),0.975,3)],clrs(i,:),0.2,'Time','$$R^{VL}(t)$$');
end
legend(hf3,{'Asx','VL','PKDL','All'})
saveas(gcf,'RtVL')
saveas(gcf,'RtVL.eps','epsc')
saveas(gcf,'RtVL.png')

%% Plot average distance to infector for VL cases and asymptomatic individuals
figure; 
PlotRjDistn(di,2,'$$\overline{d_i}$$ (m)','latex',dfltclrs(1,:))
saveas(gcf,'diDistn')
saveas(gcf,'diDistn.eps','epsc')
saveas(gcf,'diDistn.png')
median(mean(di,2,'omitnan'),'omitnan')

figure;
PlotRjDistn(di_A,2,'$$\overline{d_i}$$ (m)','latex',clrs(1,:)); hold on
PlotRjDistn(di_I,2,'$$\overline{d_i}$$ (m)','latex',clrs(2,:))
legend(gca,{'Asx','VL'})
hold off
saveas(gcf,'diDistnIandA')
saveas(gcf,'diDistnIandA.eps','epsc')
saveas(gcf,'diDistnIandA.png')
median(mean(di_A,2,'omitnan'),'omitnan')
median(mean(di_I,2,'omitnan'),'omitnan')

figure;
PlotRjDistn(diA,2,'$$\overline{d_i}$$ (m)','latex',clrs(1,:)); hold on
PlotRjDistn(diI,2,'$$\overline{d_i}$$ (m)','latex',clrs(2,:)) 
PlotRjDistn(diP,2,'$$\overline{d_i}$$ (m)','latex',clrs(3,:))
legend(gca,{'Asx','VL','PKDL'})
hold off
saveas(gcf,'diAIP')
saveas(gcf,'diAIP.eps','epsc')
saveas(gcf,'diAIP.png')

figure;
PlotRjDistn(diA_I,2,'$$\overline{d_i}$$ (m)','latex',clrs(1,:)); hold on
PlotRjDistn(diI_I,2,'$$\overline{d_i}$$ (m)','latex',clrs(2,:)) 
PlotRjDistn(diP_I,2,'$$\overline{d_i}$$ (m)','latex',clrs(3,:))
legend(gca,{'Asx','VL','PKDL'})
hold off
saveas(gcf,'di_I')
saveas(gcf,'di_I.eps','epsc')
saveas(gcf,'di_I.png')

%% Plot average distance to infector for VL cases and asymptomatic individuals
figure; 
PlotRjDistn(ti,2,'$$\overline{\tau_i}$$ (months)','latex',dfltclrs(1,:))
saveas(gcf,'tiDistn')
saveas(gcf,'tiDistn.eps','epsc')
saveas(gcf,'tiDistn.png')
min(mean(ti,2,'omitnan'))
median(mean(ti,2,'omitnan'),'omitnan')
max(mean(ti,2,'omitnan'))

figure;
PlotRjDistn(ti_A,2,'$$\overline{\tau_i}$$ (months)','latex',clrs(1,:)); hold on
PlotRjDistn(ti_I,2,'$$\overline{\tau_i}$$ (months)','latex',clrs(2,:))
legend(gca,{'Asx','VL'})
hold off
saveas(gcf,'tiDistnIandA')
saveas(gcf,'tiDistnIandA.eps','epsc')
saveas(gcf,'tiDistnIandA.png')
min(mean(ti_A,2,'omitnan'))
median(mean(ti_A,2,'omitnan'),'omitnan')
max(mean(ti_A,2,'omitnan'))
min(mean(ti_I,2,'omitnan'))
median(mean(ti_I,2,'omitnan'),'omitnan')
max(mean(ti_I,2,'omitnan'))

figure;
PlotRjDistn(tiA,2,'$$\overline{\tau_i}$$ (months)','latex',clrs(1,:)); hold on
PlotRjDistn(tiI,2,'$$\overline{\tau_i}$$ (months)','latex',clrs(2,:)) 
PlotRjDistn(tiP,2,'$$\overline{\tau_i}$$ (months)','latex',clrs(3,:))
legend(gca,{'Asx','VL','PKDL'})
hold off
saveas(gcf,'tiAIP')
saveas(gcf,'tiAIP.eps','epsc')
saveas(gcf,'tiAIP.png')

figure;
PlotRjDistn(tiA_I,2,'$$\overline{\tau_i}$$ (months)','latex',clrs(1,:)); hold on
PlotRjDistn(tiI_I,2,'$$\overline{\tau_i}$$ (months)','latex',clrs(2,:)) 
PlotRjDistn(tiP_I,2,'$$\overline{\tau_i}$$ (months)','latex',clrs(3,:))
legend(gca,{'Asx','VL','PKDL'})
hold off
saveas(gcf,'ti_I')
saveas(gcf,'ti_I.eps','epsc')
saveas(gcf,'ti_I.png')

%% Plot average distance and time to infectees for VL and PKDL infectors
srcidx=[4,6];
[meandists,meantimes]=RunPlotDistsAndTimes(infctr,src,dists,times,srcidx,clrs(2:3,:),'MeanInfctrToInfcteeDists','MeanInfctrOnsetToInfcteeInfctnTimes');
% Plot average distance and time to infectees for most likely VL and PKDL
% infectors
[meandistsmax,meantimesmax]=RunPlotDistsAndTimes(infctrmax,srcmax,distsmax,timesmax,srcidx,clrs(2:3,:),'MeanInfctrToInfcteeDistsMostLikelyInfctrs','MeanInfctrOnsetToInfcteeInfctnTimesMostLikelyInfctrs');

%% Save output
save('ReffCalcs','infctn','infctr','src','dists','times','SI','onsetinfctr','rcvryinfctr','meandists','meantimes','infctrmax','srcmax','distsmax','timesmax','SImax','onsetinfctrmax','rcvryinfctrmax','meandistsmax','meantimesmax','RjA_I','RjI_I','RjP_I','Rj_I','Rts_I','Rt_I','diA_I','diI_I','diP_I','di_I','tiA_I','tiI_I','tiP_I','ti_I','infctnA','infctrA','srcA','distsA','timesA','SIA','onsetinfctrA','rcvryinfctrA','infctrmaxA','srcmaxA','distsmaxA','timesmaxA','SImaxA','onsetinfctrmaxA','rcvryinfctrmaxA','RjA_A','RjI_A','RjP_A','Rj_A','Rts_A','Rt_A','diA_A','diI_A','diP_A','di_A','tiA_A','tiI_A','tiP_A','ti_A','iters','OT','Rj','Rts','Rt')
