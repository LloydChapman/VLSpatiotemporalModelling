function VLStmMCMC2010FastLklhdMgrtn(data,r1,p10,a,b,p2,u,beta0,alpha0,epsilon0,delta0,lambda0,h0,h1,h2,h3,hmssng,h40,pI0,typ,inclLST,niters,plotOutpt,rslts,para)
 
if nargin==0
    %%
    load('~/Dropbox/Visceral Leishmaniasis/CarynBernData/2010data/data_final2.mat') % database
    r1=3; %2; % r parameter of NegBin(r1,p1) distribution for incubation period
    mu=5; % mean incubation period in months
    p10=r1/(mu-1+r1);
    % Shape parameters for beta prior for p1 parameter of NegBin(r1,p1) distn
    a=22; %20; % 
    b=28; %38; %
    p2=1/5;
    u=1:4; %1:3; %[1:3,7];%[1:3,6]; %2:3; %[1:3,5]; %[1,3,4]; %1:2; % transmission parameters to update
    beta0=3; %1; %0.2; %6; %2; % spatial transmission rate constant
    alpha0=100; %50; % distance scale factor for spatial kernel;
    epsilon0=1e-3; %1e-4; %1e-2; % background transmission rate
    delta0=1e-3; %0; % additional within-HH transmission rate
    s1=load('~/Dropbox/Visceral Leishmaniasis/CarynBernData/2004data/SpatiotemporalModelling/data_final.mat');
    [pars,~]=FitCatModLST3(s1.data,p2);
    lambda0=pars/12;
    h0=0.03; %1/40; % relative infectiousness of pre-symptomatics
    h1=9/26/(10/15); % relative infectiousness of macular/papular PKDL cases
    h2=((9/26+18/21)/2)/(10/15); % assumed relative infectiousness of plaque PKDL cases (=(h1+h3)/2)
    h3=18/21/(10/15); % relative infectiousness of nodular PKDL cases
%     h4=0.03; %0; %0.1; % relative infectiousness of asymptomatic individuals
    h40=0.01; %0; % relative infectiousness of asymptomatic individuals
    pI0=0.15; %0.25; % proportion of infections that lead to VL
    hmssng=(101/138*h1+31/138*h2+6/138*h3); % use average relative infectiousness of PKDL cases for cases that weren't physically examined
    typ='Exp'; % type of transmission kernel
    inclLST=false; % include LST data: false=don't include, true=include
    niters=1e4; %1e5; %1e3; %5e3; %5e2; %10; %1e2; % no. of MCMC steps
    plotOutpt=true; % plot output or not (if=false)
    rslts='MCMC_NBIP_PKDL_ASX147.mat'; % name of results file for output
    para=1; %1:19; %[1:6,16:19]; %1:8; % % all paras
end
 
%% LOAD DATA 
% Select data for para
data=data(ismember(data.PARA,para),:);
 
% % Restrict dataset to one observation per individual
% [~,ic,~]=unique(data.ORIG_ID);
% data=data(ic,:);
 
% Rename longitude and latitude variables
data.Properties.VariableNames{'HHNEWLNG'}='longitude';
data.Properties.VariableNames{'HHNEWLAT'}='latitude';
 
% Number of individuals
n=size(data,1);
 
% % Set mean incubation period and treatment duration
% mu=5; % months
% FOR NOW assume cases stop being infectious shortly after commencing treatment
durRX=0; %1; % month
 
% Set start year as 2001 and start month as Jan so infection times of those
% that had KA onset in Jan 2002 can be included
startyr=2002; %2000; %2001; %1998; %
startmo=1;
 
% Set end year and month
endyr=2010;
endmo=12;
 
% Set maximum incubation period for cases with onset at start of study
maxIP=12; %6; %11; %

% Set maximum incubation period for cases with onset shortly after
% migration in
maxIP_IM=4;

% % Probability of infection leading to VL
% pI=0.15; %0 % N.B. keep fixed for now
% 
% % Probability parameter for asymptomatic infection duration distribution
% p2=1/5; % N.B. keep fixed for now
% 
% % Fixed duration of asymptomatic infection
% durA=3;

%% MAKE VECTORS OF EVENT TIMES
% Make STATA month origin
origin=stata_month(startyr,startmo)-1;

INTMIG_OUT=(data.INTMIG_OUT==1);
INTMIG_IN=(data.INTMIG_IN==1);

% Create logical index vector for duplicate observations of KA for KA cases 
% who migrated internally
KothrObs=((data.KA>=data.MIG_OUT&INTMIG_OUT)|(data.KA<data.MIG_IN&INTMIG_IN));

% Create logical index vector for duplicate observations of PKDL for PKDL 
% cases who migrated internally
PothrObs=((data.PKDL>=data.MIG_OUT&INTMIG_OUT)|(data.PKDL<data.MIG_IN&INTMIG_IN));

% Create indicator variable for KA onset between Jan 2002 and Dec 2010
K02_10=(data.KA_1210&data.KAYR>=startyr&~KothrObs);

% Create logical index vector for active KA at start of study
actvK=((data.KA<=origin&data.KARX>origin)|(data.KAYR==startyr-1&isnan(data.KARX)))&isnan(data.MIG_IN);

% Create logical index vector for pre-2002 KA
prevK=(data.KAYR<startyr&~actvK);

% Create logical index vector for KA onset before or at migration in
IpreEXTIM=(K02_10&data.EXTMIG_IN==1&data.KA<=data.MIG_IN&data.KA>origin+maxIP); 

% Create logical index vector for KA onset within 6 months of migration in
EXTIMsoonI=(K02_10&data.EXTMIG_IN==1&data.KA>=data.MIG_IN&data.KA-data.MIG_IN<=maxIP_IM&data.KA>origin+maxIP);

% Create logical index vector for KA onset before or at internal migration in
IpreINTIM=(data.KA<=data.MIG_IN&INTMIG_IN&~prevK);

% Create logical index vector for PKDL onset before or at internal migration in
PpreINTIM=(data.PKDL<=data.MIG_IN&INTMIG_IN);

% Create logical index vector for PKDL onset (without prior KA) before or at external migration in
PpreEXTIM=(isnan(data.KA)&data.PKDL<=data.MIG_IN&data.EXTMIG_IN==1);

% Create a logical index vector for KA cases that suffered treatment
% failure (include CHMP78102 who had PKDL onset 1 month after KA treatment)
RXF=(data.ALLRXFAIL==1&(data.MOS_RX_NEW_SX==0|(data.PKDL-data.KARX==1))&~KothrObs);

% Create a logical index vector for KA cases that suffered relapse (exclude
% PJNA58204 who suffered relapse after migrating out of study area)
REL=(data.ALLRXFAIL==1&~(data.MOS_RX_NEW_SX==0|(data.PKDL-data.KARX==1))&~KothrObs&~(data.KARX+data.MOS_RX_NEW_SX>data.MIG_OUT));

% Creat logical index vector for PKDL cases that were treated
RXP=(data.RXD_PKDL==1);

% Make vectors of KA onset, KA recovery, KA relapse, KA relapse recovery, 
% PKDL onset, PKDL resolution, birth, death, immigration and emigration times
tI=data.KA-origin;
% % tI(~K02_10&~prevK&~actvK)=Inf;
% % % Remove onset times for pre-2002 KA cases
% % tI(prevK)=NaN;
tR=data.KARX+durRX-origin;
% Add 1 month of infectiousness for treatment failure KA cases
tR(RXF)=tR(RXF)+1;
% % tR(~K02_10&~prevK&~actvK)=Inf;
% % % Remove treatment times for pre-2002 KA cases
% % tR(prevK)=NaN;
% tRL=NaN(n,1);
tRL=data.KARX+data.MOS_RX_NEW_SX-origin;
tRL(RXF)=NaN;
tRLR=NaN(n,1);
tP=data.PKDL-origin;
tRP=data.PKDL+data.PKDL_DUR-origin;
% Overwrite resolution times with treatment times for treated PKDL cases.
tRP(RXP)=max(data.PKRX(RXP),data.PKRX2(RXP))-origin; % use 2nd treatment time for case with 2 PKDL treatments
tB=data.DOB-origin;
% Give individuals missing DOB dummy DOB of 0
% tB(isnan(data.DOB))=0;
tD=data.DEATH-origin;
tIM=data.MIG_IN-origin;
tEM=data.MIG_OUT-origin;
% tEN=data.ENTRY-origin;
% tEX=data.EXIT-origin;
% tD(isnan(data.DEATH_Y))=Inf;
% End time
tmax=stata_month(endyr,endmo)-origin;
 
% Make index vectors for KA cases with and without onset times (O and NO), 
% treated indvdls, cases missing treatment times (NR), LST positives, 
% relapsers, treated relapsers, and births and deaths
OR=find(K02_10&~isnan(tI)&(~isnan(tR)|~isnan(tD))); % KA: onset and treatment/death time
NONR=find(K02_10&isnan(tI)&isnan(tR)); % KA: no onset or treatment time
ONR=find(K02_10&~isnan(tI)&isnan(tR)&isnan(tD)); % KA: onset but no treatment time & didn't die from KA
RNO=find(K02_10&isnan(tI)&~isnan(tR)); % KA: treatment but no onset time
NO=[NONR;RNO];
I=sort([OR;NONR;ONR;RNO]);
if inclLST
    NI=(~prevK&~K02_10);
    L02=find(data.LST02==1&NI&isnan(tB)); % LST+ in 2002, excluding LST+ves who have had KA, 1 individual who has KA and is LST+, and 1 individual who gets KA later
    L02B=find(data.LST02==1&NI&~isnan(tB)); % LST+ in 2002, born after START time
    L03C=find(data.LST03==1&data.LST02==0&NI); % LST conversion between 2002 and 2003
    L03=find(data.LST03==1&isnan(data.LST02)&NI&isnan(tB)); % LST+ in 2003, no LST 2002
    L03B=find(data.LST03==1&isnan(data.LST02)&NI&~isnan(tB)); % LST+ in 2003, no LST in 2002, born after START time
    L04C=find(data.LST04==1&data.LST02==0&data.LST03==0&NI); % LST conversion between 2003 and 2004
    % N.B. No pstve LST04 asymptomatic individuals with missing LST02 or LST03
    L=sort([L02;L02B;L03C;L03;L03B;L04C]);
end
% RL=find(~isnan(tRL));
% RLR=find(~isnan(tRLR));
RLO=find(REL&~isnan(tRL));
RLNO=find(REL&isnan(tRL));
RL=[RLO;RLNO];
P=find(~isnan(tP)&~PothrObs);
% RP=find(~isnan(tRP));
IandP=sort(intersect(I,P));
% INP=sort(setdiff(I,P));
PIA=setdiff(setdiff(P,I),find(PpreEXTIM)); % PKDL w/o KA during study
PI=intersect(PIA,union(find(actvK|prevK),find(tP>tIM&tI<tIM))); % PKDL with KA onset before startyr or before immigration
PA=setdiff(PIA,PI); % PKDL without prior KA (N.B. none of these PKDL cases internally migrated, so don't need to worry about internal migration with asymptomatic infection/PKDL for them)
% PI=setdiff(P,I); % PKDL w/o KA during study
B=find(tB>0); % exclude people estimated to have been born before START time
D=find(~isnan(tD));
DpreR=find(~prevK&~isnan(tI)&isnan(tR)&~isnan(tD));
RpreD=setdiff(I,DpreR);
IM=find(~isnan(tIM));
EM=find(~isnan(tEM));
IMI=find(KothrObs&tI<tIM&tR>tIM);
IMP=find(PothrObs&tP<tIM&tRP>tIM|PpreEXTIM);
% Create index vectors for 1st and 2nd observations of internal migrators 
% (complicated definition to allow for subsetting of the data and to ensure
% 1st and 2nd obs are matched in IM_OUT and IM_IN)
IM_OUT=find(ismember(data.RESP_ID,data.ORIG_ID(INTMIG_IN))&INTMIG_OUT); 
IM_IN=find(ismember(data.ORIG_ID,data.RESP_ID(INTMIG_OUT))&INTMIG_IN);
% EN=find(~isnan(tEN));
% EX=find(~isnan(tEX));

AOR=find(actvK&~isnan(tI)&~isnan(tR));
AONR=find(actvK&~isnan(tI)&isnan(tR));
ANONR=find(actvK&isnan(tI)&isnan(tR));
A=sort([AOR;AONR;ANONR]);
IPNIA=[I;PI;A;IMI;IMP];

% Number of KA cases
nI=numel(I);
% Number of cases with both onset and treatment times
nOR=numel(OR);
% Number of KA cases w/o onset or treatment times
nNONR=numel(NONR);
% Number of KA cases w/ onset time but no treatment time
nONR=numel(ONR);
% Number of KA cases w/ treatment time but no onset time
nRNO=numel(RNO);
% Number of KA cases w/ no onset time
nNO=numel(NO);
% Number of relapse KA cases w/ second onset time
nRL=numel(RL);
% Number of relapse KA cases w/ second onset time
nRLO=numel(RLO);
% Number of relapse KA cases w/o second onset time
nRLNO=numel(RLNO);
% Number of KA cases w/ onset between startyr and endyr and later PKDL
nIandP=numel(IandP);
% Number of internal immigrants with KA/PKDL at time of migration
nIMI=numel(IMI);
nIMP=numel(IMP);
% % Number of KA cases w/o later PKDL
% nINP=numel(INP);
% Number of PKDL cases w/ KA onset before startyr or before immigration
nPI=numel(PI);
% Number of PKDL cases w/o prior KA
nPA=numel(PA);
% Number of cases who may have had active KA at start of study
nA=numel(A);
% Number of potentially active KA cases without treatment time
nAONR=numel(AONR);
% Number of potentially active KA cases without onset or treatment time
nANONR=numel(ANONR);
% Number of potential infection sources
nIPNIA=numel(IPNIA);
% Number of individuals born in first month of study
nB1=sum(tB==1);

% % Probability of asymptomatic infection progressing to dormant infection
% pD=16/round((1-pI)/pI*nI);

%% DRAW INITIAL MISSING ONSET AND TREATMENT TIMES
% Create vectors of lower and upper bounds for onset month
tIlb=NaN(n,1);
tIlb(I)=max(stata_month(data.KAYR(I),1)-origin,tB(I)+2); % +2 because individuals can only be infected 1 month after birth or immigration
tIub=NaN(n,1);
% ADD tP AND tR IN HERE AND REMOVE FORM CONDTNS BELOW!
tIub(I)=min(min(min(stata_month(data.KAYR(I),12)-origin,tR(I)-1),tP(I)-2),tD(I)-1);
% Fit negative binomial distribution to onset-to-treatment (OT) times
[r0,p0]=FitOTdistn(tI,tR); % use all cases with onset and treatment times (not only those with onset in 2002-2010) as this distn is mostly used to impute missing times for active KA cases at start of study

% For individuals with missing onset, diagnosis and treatment times draw
% onset time at random from onset year
for i=1:nNONR
%     while isnan(tI(NONR(i))) || tI(NONR(i))>=tD(NONR(i))
%         tI(NONR(i))=(data.KAYR(NONR(i))-startyr)*12+randi(12)-startmo+1;
%     end
    tI(NONR(i))=randi([tIlb(NONR(i)),tIub(NONR(i))],1);
end
% Draw treatment times for these individuals using OT distn
for i=1:nNONR
    while isnan(tR(NONR(i))) || tR(NONR(i))>min(tD(NONR(i)),tmax)
        tR(NONR(i))=tI(NONR(i))+nbinrnd(r0,p0)+1;
    end
end
% For individuals with onset time but no treatment time draw treatment time
% using OT distn
for i=1:nONR
    while isnan(tR(ONR(i))) || tR(ONR(i))>min(tD(ONR(i)),tmax)
        tR(ONR(i))=tI(ONR(i))+nbinrnd(r0,p0)+1;
    end
end
% For individuals with treatment time but no onset time draw onset time at
% random from onset year
for i=1:nRNO
%     while isnan(tI(RNO(i))) || tI(RNO(i))>=tR(RNO(i))
%         tI(RNO(i))=(data.KAYR(RNO(i))-startyr)*12+randi(12)-startmo+1;
%     end
    tI(RNO(i))=randi([tIlb(RNO(i)),tIub(RNO(i))],1);
end

% Create vectors of lower and upper bounds for onset month for active KA
% cases at start of study
tIlbA=NaN(n,1);
tIlbA(A)=max(stata_month(data.KAYR(A),1)-origin,tB(A)+2);
tIubA=NaN(n,1);
tIubA(A)=min(min(min(stata_month(data.KAYR(A),12)-origin,tR(A)-1),tP(A)-2),tD(A)-1);

for i=1:nANONR
%     while isnan(tI(ANONR(i))) || tI(ANONR(i))>=tD(ANONR(i))
%         tI(ANONR(i))=(data.KAYR(ANONR(i))-startyr)*12+randi(12)-startmo+1;
%     end
    tI(ANONR(i))=randi([tIlbA(ANONR(i)),tIubA(ANONR(i))],1);
end
for i=1:nANONR
    while isnan(tR(ANONR(i))) || tR(ANONR(i))>min(tD(ANONR(i)),tmax)
        tR(ANONR(i))=tI(ANONR(i))+nbinrnd(r0,p0)+1;
    end
end
for i=1:nAONR
    while isnan(tR(AONR(i))) || tR(AONR(i))>min(tD(AONR(i)),tmax)
        tR(AONR(i))=tI(AONR(i))+nbinrnd(r0,p0)+1;
    end
end

% Make vector of recovery times for treated individuals and death times for
% untreated individuals
tRorD=tR;
tRorD(DpreR)=tD(DpreR);

%% DRAW INITIAL RELAPSE AND RELAPSE TREATMENT TIMES
p4=mle(tRL(RLO)-tR(RLO)-1,'distribution','geo');
% x=1:30; figure; histogram(tRL(RLO)-tR(RLO),3,'Norm','pdf'); hold on; plot(x,geopdf(x-1,p4)); hold off

for i=1:nRLNO
    j=RLNO(i);
    while isnan(tRL(j)) || tRL(j)>min(min(min(tEM(j),tP(j)),tD(j))-1,tmax)-1
        tRL(j)=tR(j)+geornd(p4)+1;
    end
end

for i=1:nRL
    j=RL(i);
    while isnan(tRLR(j)) || tRLR(j)>min(min(min(tEM(j),tP(j)),tD(j))-1,tmax)
        tRLR(j)=tRL(j)+nbinrnd(r0,p0)+1;
    end
end

%% DRAW INITIAL PRE-SYMPTOMATIC INFECTION TIMES
% Draw infection times from negative binomial distribution with parameters
% r1 and p1
tE=Inf(n,1); % initialise infection time vector
% tE(prevK)=0; % set dummy infection time of 0 for previous KA cases
tE(I)=-Inf; % set infection times so that they will all be updated below
% % Make vector of birth times for KA cases born during study
% tBI=NaN(n,1);
% tBI(setdiff(I,B))=0; % for those not born during study set dummy birth time of 0
% tBI(intersect(I,B))=tB(intersect(I,B));
for i=1:nI
    % loop until infection time is not before month after birth and cases 
    % with onset after initial window have infection times after first month
    j=I(i);
    while tE(j)<tB(j)+1 || (tI(j)>maxIP && tE(j)<1) || (tI(j)-tIM(j)-1>maxIP_IM && tE(j)<tIM(j)+1) % individuals can't be infected until month after birth/immigration
        tE(j)=tI(j)-(nbinrnd(r1,p10)+1);
    end
end
% figure; histogram(tE(I),'BinMethod','Integers')

% Make index vector to exclude cases with onset before maxIP or before or shortly after migration in
I1=find(tI>maxIP&~IpreEXTIM&~EXTIMsoonI&~KothrObs&~isinf(tI)); % ~isinf(tI) is redundant with redefinition of tI

%% INFECTIOUS PRESSURE FROM KA AND PKDL CASES
% Calculate distances between HHs
[dHH,ia,ib]=CalcHHDists(data);
nHH=numel(ia);
f=histcounts(ib,1:nHH+1)';
if strcmp(typ,'Cauchy')
    dHHsqrd=dHH.^2;
else
    dHHsqrd=[];
end

[KHH,K0old]=Knl_fast(dHH,dHHsqrd,alpha0,beta0,typ,n,nHH,ib,f,1:n);

% d0=bsxfun(@eq,(1:nHH)',ib');
% d0=sparse(d0);
d0=speye(nHH);
rateHHA=beta0*KHH+delta0*d0;
rateHH=rateHHA(:,ib(IPNIA));

% Infectiousness over time
% Make vector of lesion-specific infectiousnesses for PKDL cases
hv=NaN(n,1);
hv(strcmp(data.TYPE,'MAC_PAP'))=h1;
hv(strcmp(data.TYPE,'PLQ'))=h2;
hv(strcmp(data.TYPE,'NOD'))=h3;
hv(~isnan(tP)&strcmp(data.TYPE,''))=hmssng;

h=zeros(nIPNIA,tmax); % initialise infectiousness matrix
for i=1:nIandP
    j=IandP(i);
    h(I==j,max(tP(j),tIM(j))+1:min(min(tRP(j),tEM(j)),tmax))=hv(j);
end
for i=1:nI
    j=I(i);
    h(i,max(max(0,tIM(j)),tE(j))+1:tI(j))=h0; % infectiousness from month of infection up to month before symptom onset
    h(i,max(tIM(j),tI(j))+1:min(tEM(j),tRorD(j)))=1; % infectiousness up to month before month of treatment or death
end
for i=1:nPI
    j=PI(i);
    h(nI+i,tP(j)+1:min(min(tRP(j),tEM(j)),tmax))=hv(j);
end
% Assume individuals return to full infectiousness on relapse until
% re-treatment
for i=1:nRL
    j=RL(i);
    h(IPNIA==j,tRL(j)+1:min(min(tRLR(j),tEM(j)),tmax))=1;
end
for i=1:nA
    j=A(i);
    h(nI+nPI+i,1:tRorD(j))=1;
end
for i=1:nIMI
    j=IMI(i);
    h(nI+nPI+nA+i,tIM(j)+1:tRorD(j))=1;
%     h(j,tIM(j):tRorD(j)-1)=1;
end
for i=1:nIMP
    j=IMP(i);
    h(nI+nPI+nA+nIMI+i,tIM(j)+1:min(min(tRP(j),tEM(j)),tmax))=hv(j);
%     h(j,tIM(j):min(tRP(j)-1,tmax))=hv(j);
end
h=sparse(h);

% lambdaHHI=rateHH*h+epsilon0; % HH-level infectious pressure
lambdaHHI=rateHH*h; % HH-level infectious pressure
lambdaI=lambdaHHI(ib,:); % expand to individual-level infectious pressure

%% DRAW INITIAL ASYMPTOMATIC INFECTION AND RECOVERY TIMES
% SusA=setdiff((1:n)',union(IPNIA,find(actvK|prevK|IpreEXTIM|EXTIMsoonI|IpreINTIM|PpreINTIM|PpreEXTIM|KothrObs|PothrObs)));
SusA=setdiff((1:n)',union(union(IPNIA,PA),find(actvK|prevK|IpreEXTIM|EXTIMsoonI|IpreINTIM|PpreINTIM|PpreEXTIM|KothrObs|PothrObs)));
nSusA=numel(SusA);
% s1=load('~/Dropbox/Visceral Leishmaniasis/CarynBernData/2004data/SpatiotemporalModelling/data_final.mat');
% [pars,~]=FitCatModLST3(s1.data,p2);
% % age=max((1-tB)/12,0); % Should this be max(-tB/12,0)?
age=max(-tB,0); 
rng=[max(max(tB,tIM)+1,0) min(min(tEM,tD),tmax+1)];
% Index vector for individuals present at (born or imigrated before) time 0
Pres0=find(rng(:,1)==0);
% figure; histogram(rng(:,1),'BinM','int'); hold on; histogram(rng(:,2),'BinM','int'); hold off
SusA0=intersect(SusA,Pres0);
% lambda0=pars/12; %0.001/12; %0.01/12; %0.001/12;
prob0=ProbInitStatus(age,lambda0,p2);
% % prob0=ProbInitStatus(age,0.001/12,pI,p2);
prob0=[prob0,1-sum(prob0,2)];
prob0(rng(:,1)>0,1)=1;
prob0(rng(:,1)>0,2:3)=0;
Stat0=NaN(n,1);
for i=1:numel(SusA0)
    j=SusA0(i);
    Stat0(j)=randsample(1:3,1,true,prob0(j,:));
end
% probLpos=PrevVA(age,pars(1),pars(2));
% prevA=SusA0(rand(numel(SusA0),1)<probLpos(SusA0));%0.3); %
% Sus=setdiff(SusA,prevA);
actvA=find(Stat0==2);
nactvA=numel(actvA);
prevA=find(Stat0==3);
IM_INprevAactvA=IM_IN(ismember(IM_OUT,[prevA;actvA]));
Sus=setdiff(SusA,[union(prevA,actvA);IM_INprevAactvA]);
% For internal migrators, randomly pick which observation to remove
IM_OUT_IN=[IM_OUT,IM_IN];
IM_OUTSus=intersect(IM_OUT,Sus);
nIM_OUTSus=numel(IM_OUTSus);
% rmIM=IM_OUT_IN(sub2ind(size(IM_OUT_IN),find(ismember(IM_OUT,IM_OUTSus)),randi(2,numel(IM_OUTSus),1)));
% kpIM=setdiff(IM_OUT_IN(ismember(IM_OUT,IM_OUTSus),:),rmIM);
rmIM=NaN(nIM_OUTSus,1);
kpIM=NaN(nIM_OUTSus,1);
for i=1:nIM_OUTSus
    rmIM(i)=IM_OUT_IN(IM_OUT==IM_OUTSus(i),randi(2));
    kpIM(i)=IM_OUT_IN(IM_OUT==IM_OUTSus(i),IM_OUT_IN(IM_OUT==IM_OUTSus(i),:)~=rmIM(i));
end
Sus=setdiff(Sus,rmIM);
nSus=numel(Sus);
% Calculate cumulative infectious pressure
% LambdaHH=sum(lambdaHH,2);
% Lambda=LambdaHH(ib,:);
rngm=ones(n,tmax);
for i=1:n
    if rng(i,1)>1 || rng(i,2)<tmax+1
        rngm(i,[1:rng(i,1)-1,rng(i,2):tmax])=0;
    end
end
% figure; spy(1-rngm(1:100,:)); axis square
% LambdaI=sum(lambdaI.*rngm,2);
% CumLambda=cumsum(lambda.*rngm,2);
% CumLambdaI=cumsum(lambdaI,2);
% CumLambdaI=[zeros(n,1),CumLambdaI(:,1:end-1)];
CumLambdaI=cumsum(lambdaI.*rngm,2);
CumLambdaI=[zeros(n,1),CumLambdaI(:,1:end-1)];
% CumLambda=[zeros(n,1),CumLambda];
% probA=exp(-(1-pI)*CumLambda).*(1-exp(-(1-pI)*lambda));
% probA=[1-prob0(:,1),bsxfun(@times,prob0(:,1),[exp(-(1-pI)*CumLambdaI).*(1-exp(-(1-pI)*lambdaI)).*rngm,1-sum(exp(-(1-pI)*CumLambdaI).*(1-exp(-(1-pI)*lambdaI)).*rngm,2)])]; %PRETTY SURE THIS IS WRONG
% probA=[1-prob0(:,1),bsxfun(@times,prob0(:,1),[exp(-CumLambdaI).*(1-exp(-(1-pI)*lambdaI)).*rngm,1-sum(exp(-CumLambdaI).*(1-exp(-(1-pI)*lambdaI)).*rngm,2)])];
probA1=[exp(-CumLambdaI).*(1-exp(-(1-pI0)*lambdaI.*rngm)),exp(-sum(lambdaI.*rngm,2))];
probA1=bsxfun(@rdivide,probA1,sum(probA1,2));
probA=[1-prob0(:,1),bsxfun(@times,prob0(:,1),probA1)];
% % probA=[40*(rng(:,1)==0),rngm,80*(rng(:,2)==tmax+1)];
% % probA=[40*(rng(:,1)==0),rngm,80*ones(n,1)];
% % probA=[40*(rng(:,1)==0),rngm,40*ones(n,1)];
% probA=[80*(rng(:,1)==0),rngm,80*ones(n,1)];
% % cprobA=bsxfun(@rdivide,cumsum(probA,2),sum(probA,2));
% % 
LambdaI=sum(probA(:,1:end-1),2);


% % probA=[rng(:,1)==0,rngm,ones(n,1)];
% probA=bsxfun(@rdivide,probA,sum(probA,2));

% probA=[100*(1-prob0(:,1)),rngm,80*(rng(:,2)==tmax+1)];
% probA=[1-prob0(:,1),rngm/100,prob0(:,1).*(1-sum(exp(-CumLambdaI).*(1-exp(-(1-pI)*lambdaI)).*rngm,2))];
% probA=[15*(rng(:,1)==0),rngm,80*(rng(:,2)==tmax+1)];
% probA(PA,end)=0;
% probA(PA,:)=bsxfun(@rdivide,probA(PA,:),sum(probA(PA,:),2));

% Pick asymptomatic individuals with probability proportional to the
% cumulative infectious pressure on them during the study
% nAsx=round((1-pI)/pI*nI);
tA=NaN(n,1);
tRA=NaN(n,1);
tA(prevA)=0;
tRA(prevA)=0;
tA(actvA)=0;
for i=1:numel(actvA)
    j=actvA(i);
    if ismember(j,IM_OUT)
        j1=IM_IN(IM_OUT==j);
        probAIP=[geopdf(0:rng(j1,2)-2,p2),1-geocdf(rng(j1,2)-2,p2)];
        fwdRA=[1:rng(j1,2)-1,tmax+1];
    else
        probAIP=[geopdf(0:rng(j,2)-2,p2),1-geocdf(rng(j,2)-2,p2)];
        fwdRA=[1:rng(j,2)-1,tmax+1];
    end    
    tRA(j)=fwdRA(randsample(numel(fwdRA),1,true,probAIP));
end
tA(IM_INprevAactvA)=tmax+2;
tRA(IM_INprevAactvA)=tmax+2;
IM_OUTactvA=IM_OUT(ismember(IM_OUT,actvA));
IM_INactvA=IM_IN(ismember(IM_OUT,actvA));
RAobs1actvA=IM_OUTactvA(tRA(IM_OUTactvA)>rng(IM_OUTactvA,2)-1);
RAobs2actvA=IM_INactvA(tRA(IM_OUTactvA)>rng(IM_OUTactvA,2)-1);
% for i=1:nSus
%     j=Sus(i);
%     tA(j)=randsample(1:tmax+1,1,true,probA(j,2:tmax+2));
%     probAIP=[geopdf(0:rng(j,2)-tA(j)-2,p2),1-geocdf(rng(j,2)-tA(j)-2,p2)];
%     fwdRA=[tA(j)+1:rng(j,2)-1,tmax+1];
%     tRA(j)=fwdRA(randsample(numel(fwdRA),1,true,probAIP));
% end
% Asx=find(tA>0 & tA<tmax+1);

% nAsx=round((1-pI)/pI*nI)-nPA; %0;
nAsx=min(round((1-pI0)/pI0*nI)-nPA,nSus);
Asx=Sus(datasample(1:nSus,nAsx,'Replace',false,'Weights',LambdaI(Sus))); %)); %
% Asx=Sus(sort(randperm(nSus,nAsx)));
% Make index vector of individuals still susceptible at end of study
Susend=[setdiff(Sus,Asx);rmIM(~ismember(kpIM,Asx));IM_OUT(ismember(IM_IN,Asx))];
tA(Susend)=tmax+1;
tRA(Susend)=tmax+1;
% Set asymptomatic infection and recovery time of 2nd observation of
% internal migrators asymptomatically infected during 1st observation to
% dummy time of tmax+2
tA(IM_IN(ismember(IM_OUT,Asx)))=tmax+2;
tRA(IM_IN(ismember(IM_OUT,Asx)))=tmax+2;
for i=1:nAsx
    j=Asx(i);
%     tA(j)=randsample(max(1,rng(j,1)):rng(j,2),1,true,lambda(j,max(1,rng(j,1)):rng(j,2))); %);%
%     tA(j)=randsample(max(1,rng(j,1)):rng(j,2)-1,1,true,probA(j,(max(1,rng(j,1)):rng(j,2)-1)+1)); %);%
    tA(j)=randsample(1:tmax,1,true,probA(j,2:tmax+1));
    if ismember(j,IM_OUT)
        j1=IM_IN(IM_OUT==j);
        probAIP=[geopdf(0:rng(j1,2)-tA(j)-2,p2),1-geocdf(rng(j1,2)-tA(j)-2,p2)];
        fwdRA=[tA(j)+1:rng(j1,2)-1,tmax+1];
    else
        probAIP=[geopdf(0:rng(j,2)-tA(j)-2,p2),1-geocdf(rng(j,2)-tA(j)-2,p2)];
        fwdRA=[tA(j)+1:rng(j,2)-1,tmax+1];
    end
    tRA(j)=fwdRA(randsample(numel(fwdRA),1,true,probAIP));
end
IM_OUTAsx=IM_OUT(ismember(IM_OUT,Asx));
IM_INAsx=IM_IN(ismember(IM_OUT,Asx));
RAobs1=IM_OUTAsx(tRA(IM_OUTAsx)>rng(IM_OUTAsx,2)-1);
RAobs2=IM_INAsx(tRA(IM_OUTAsx)>rng(IM_OUTAsx,2)-1);

% save('tA','tA','tRA','prevA','actvA','Asx','Sus')
% load('tA')

% % Start chain with all non-symptomatic individuals susceptible
% tA(~isnan(tA))=tmax+1;
% tRA(~isnan(tA))=tmax+1;
% actvA=[];
% prevA=[];
% Asx=[];

RP=tP-tR;
pars1=nbinfit(RP(RP>=0));
r3=pars1(1);
p3=pars1(2);
% N.B. None of the PKDL cases without prior KA internally migrated, so
% don't need to worry about asymptomatic infection overlapping separate
% observations
for i=1:nPA
    j=PA(i);
    while isnan(tRA(j)) || tA(j)<max(1,tB(j)+1) %tA(j)<tB(j) %tRA(j)<tB(j)+1
        tRA(j)=tP(j)-nbinrnd(r3,p3);
%     end
%     if tRA(j)<=0
%         tRA(j)=0;
%         tA(j)=0;
%     else
        tA(j)=max(0,tRA(j)-(geornd(p2)+1));
    end
end
S0PA=PA(tA(PA)>0);
actvAPA=PA(tA(PA)<=0 & tRA(PA)>0);
prevAPA=PA(tRA(PA)<=0);

% figure; histogram(tA(tA>0 & tA<tmax+1),'BinMethod','Integers')
% hold on; histogram(tRA(tRA>0 & tRA<tmax+1),'BinMethod','Integers'); hold off
% figure; histogram(tRA(PA)-tA(PA),'BinM','int')

A1=[Asx;actvA];
nA1=numel(A1);

%% INFECTIOUSNESS
% Calculate initial incubation periods
IPs=zeros(nI,niters); % matrix of infection periods
IPold=NaN(n,1);
IPold(I)=tI(I)-tE(I); % SHOULD I EXCLUDE INDIVIDUALS WITH ONSET BEFORE maxIP FROM THE INCUBATION PERIOD PART OF THE LIKELIHOOD? PROBABLY
 
% % Infectiousness over time
% % Make vector of lesion-specific infectiousnesses for PKDL cases
% hv=NaN(n,1);
% hv(strcmp(data.TYPE,'MAC_PAP'))=h1;
% hv(strcmp(data.TYPE,'PLQ'))=h2;
% hv(strcmp(data.TYPE,'NOD'))=h3;
% hv(~isnan(tP)&strcmp(data.TYPE,''))=hmssng;
% % figure; plot(find(~isnan(hv))',hv(~isnan(hv)),'x')
% 
% h=zeros(nIPNIA,tmax); % initialise infectiousness matrix
% hA=zeros(nA1,tmax); % initalise asx infectiousness matrix
% h=zeros(n,tmax); % size of infectiousness matrix
% Fill in PKDL infectiousness first so that 1 case who had simultaneous
% KA and PKDL has KA infectiousness until KA treatment

%%%% N.B. Column indices of h now refer to the previous time step, i.e.
%%%% column j refers to time j-1


% for i=1:nIandP
%     j=IandP(i);
%     h(ismember(I,j),max(tP(j),tIM(j))+1:min(min(tRP(j),tEM(j)),tmax))=hv(j);
% %     h(j,max(tP(j),tIM(j)):min(min(tRP(j)-1,tEM(j)-1),tmax))=hv(j);
% end
% for i=1:nI
%     j=I(i);
%     h(i,max(max(0,tIM(j)),tE(j))+1:tI(j))=h0; % infectiousness from month of infection up to month before symptom onset
%     h(i,max(tIM(j),tI(j))+1:min(tEM(j),tRorD(j)))=1; % infectiousness up to month before month of treatment or death
% %     h(j,max(max(1,tIM(j)),tE(j)):tI(j)-1)=h0; % infectiousness from month of infection up to month before symptom onset
% %     h(j,max(tIM(j),tI(j)):min(tEM(j)-1,tRorD(j)-1))=1; % infectiousness up to month before month of treatment or death
% end
% for i=1:nPNI
%     j=PNI(i);
%     h(nI+i,tP(j)+1:min(min(tRP(j),tEM(j)),tmax))=hv(j);
% %     h(j,tP(j):min(min(tRP(j)-1,tEM(j)-1),tmax))=hv(j);
% end
% % % Assume individuals return to full infectiousness on relapse until
% % % re-treatment
% % for j=1:numel(RL)
% %     h(RL(j),tRL(RL(j)):min(tRLR(RL(j))-1,tmax))=1;
% % end
% for i=1:nA
%     j=A(i);
%     h(nI+nPNI+i,1:tRorD(j))=1;
% %     h(j,1:tRorD(j)-1)=1;
% end
% for i=1:nIMI
%     j=IMI(i);
%     h(nI+nPNI+nA+i,tIM(j)+1:tRorD(j))=1;
% %     h(j,tIM(j):tRorD(j)-1)=1;
% end
% for i=1:nIMP
%     j=IMP(i);
%     h(nI+nPNI+nA+nIMI+i,tIM(j)+1:min(tRP(j),tmax))=hv(j);
% %     h(j,tIM(j):min(tRP(j)-1,tmax))=hv(j);
% end
% for i=1:nAsx
%     j=Asx(i);
%     hA(i,tA(j)+1:min(min(min(tRA(j)-1,tEM(j)),tD(j))+1,tmax))=h4;
% end
% for i=1:nactvA
%     j=actvA(i);
%     hA(nAsx+i,1:min(min(min(tRA(j)-1,tEM(j)),tD(j))+1,tmax))=h4;
% end
% for i=1:nA1
%     j=A1(i);
%     hA(i,tA(j)+1:min(min(min(tRA(j),tEM(j)),tD(j)),tmax))=h4;
% end

hA=zeros(n,tmax);
for i=1:nA1
    j=A1(i);
    hA(j,tA(j)+1:min(min(min(tRA(j),tEM(j)),tD(j)),tmax))=h40;
end
for i=1:numel(RAobs2actvA)
    j=RAobs1actvA(i);
    j1=RAobs2actvA(i);
    hA(j1,rng(j,2)+1:min(min(tRA(j),rng(j1,2)),tmax))=h40;
end 
for i=1:numel(RAobs2)
    j=RAobs1(i);
    j1=RAobs2(i);
    hA(j1,rng(j,2)+1:min(min(tRA(j),rng(j1,2)),tmax))=h40;
end
hA=sparse(hA);
% figure; spy(hA(RAobs2,:)); axis square

hPA=zeros(nPA,tmax);
for i=1:nPA
    j=PA(i);
    hPA(i,max(1,tA(j)+1):tRA(j))=h40;
    hPA(i,tP(j)+1:min(min(tRP(j),tEM(j)),tmax))=hv(j);
end
% figure; spy(h); axis square
% figure; spy(hA); axis square
% figure; spy(hPA); axis square

%% INFECTIOUS PRESSURE
% % tstart=tic;
% % Calculate distances between HHs
% [dHH,ia,ib]=CalcHHDists(data);
% nHH=numel(ia);
% f=histcounts(ib,1:nHH+1)';
% if strcmp(typ,'Cauchy')
%     dHHsqrd=dHH.^2;
% else
%     dHHsqrd=[];
% end
 
% Calculate infectious pressure over time
% tstart=tic;
% profile on
% [KHH,K0old]=Knl_fast(dHH,dHHsqrd,alpha0,beta0,typ,n,nHH,ib,f,[IPNIA;A1]); % spatial kernel
% [KHH,K0old]=Knl_fast2(dHH,dHHsqrd,alpha0,beta0,typ,n,nHH,ib,f,IPNIA); % spatial kernel
% profile viewer
% toc(tstart)
 
% tstart=tic;
% profile on
% [K,K0old]=Knl_fast2(d,dsqrd,alpha0,typ,n); % spatial kernel
% profile viewer
% toc(tstart)
 
% tstart=tic;
% [K,K0old]=Knl(d,alpha0,typ,n); % spatial kernel
% toc(tstart)
 
% tstart=tic;
% % Calculate distances between individuals
% d=CalcDists(data);

% % Calculate infectious pressure over time
% [K,K0old]=Knl(d,alpha0,typ,n); % spatial kernel
% d0=double(d==0); % indicator matrix of individuals living in the same HHs
% d0(1:n+1:end)=0; % set diagonal (same individual) entries to 0
% d0=bsxfun(@eq,(1:nHH)',ib');
% rateHH=beta0*KHH(:,[IPNIA;A1])+delta0*d0(:,[IPNIA;A1]);
% rateHH=beta0*KHH(:,IPNIA)+delta0*d0(:,IPNIA);
% rateHHA=beta0*KHH(:,A1)+delta0*d0(:,A1);
% rateHHA=beta0*KHH+delta0*d0;
% rateHH=rateHHA(:,IPNIA);
% rate=beta0*K+delta0*d0; % calculate transmission rate, adding within-HH contribution
% toc(tstart)
 
% % NEEDS FIXING
% if ismember(4,u)
%     d0c=cell(nHH,1);
%     for i=1:nHH
%         d0c{i}=ones(f(i));
%     end
%     d0=blkdiag(d0c{:}); % indicator matrix of individuals living in the same HHs
%     d0(1:n+1:end)=0; % set diagonal (same individual) entries to 0
%     d0=d0(:,I);
%     rate=rate+delta0*d0; % calculate transmission rate, adding within-HH contribution
% end
% lambda=rate(:,I)*h(I,:)+epsilon0; % infectious pressure
% tic;
% lambdaHH=rateHH*[h;hA]+epsilon0; % HH-level infectious pressure
lambdaHHA=rateHHA(:,ib([A1;RAobs2actvA;RAobs2]))*hA([A1;RAobs2actvA;RAobs2],:);
rateHHPA=rateHHA(:,ib(PA));
lambdaHHPA=rateHHPA*hPA;
lambdaHH=rateHH*h+lambdaHHPA+lambdaHHA+epsilon0; % HH-level infectious pressure
lambda=lambdaHH(ib,:); % expand to individual-level infectious pressure
% szl=[n tmax];
% szr=[nHH nIPNIA+nAsx+nactvA];
% szh=[nIPNIA tmax];
% hidx1=find(ismember(IPNIA,I1));
% hidxA=(1:nAsx)';
% szhA=[nAsx+nactvA tmax];
% lambda(sub2ind(szl,[I1;Asx],[tE(I1);tA(Asx)]))=lambda(sub2ind(szl,[I1;Asx],[tE(I1);tA(Asx)]))-rateHH(sub2ind(szr,ib([I1;Asx]),[hidx1;nIPNIA+hidxA])).*[h(sub2ind(szh,hidx1,tE(I1)));hA(sub2ind(szhA,hidxA,tA(Asx)))]; % remove infectious pressures of cases on themselves
% lambda=repelem(lambda,f,1);
% toc

% rate=rateHH(ib,:);
% rate(sub2ind([n nIPNIA],I1,hidx1))=0;
% lambda1=rate*h+epsilon0;
% sum(sum(lambda.*S-lambda1.*S))

% d=dHH(ib,ib(IPNIA));
% d0=double(d==0);
% rate=beta0*KHH(ib,:)+delta0*d0;
% rate(sub2ind([n nIPNIA],I1,hidx1))=0;
% lambda=rate*h+epsilon0;

%% MAKE STATUS MATRICES
% Make logical matrices for pre-symptomatic infection, asymptomatic 
% infection, relapse, relapse recovery, PKDL, onset, birth, death, 
% immigration and emigration times

% Pre-symptomatic infection
tEm=false(n,tmax);
tEm((tE(I1)-1)*n+I1)=1;
tEm=sparse(tEm);
% for i=1:nI
%     j=I(i);
%     if tI(j)>maxIP && ~IpreEXTIM(j) && ~EXTIMsoonI(j) % exclude cases with onset before maxIP or before or shortly after migration in
%         tEm(j,tE(j))=1;
%     end
% end
% Asymptomatic infection
tAm=false(n,tmax);
% tAm((tA(Asx)-1)*n+Asx)=1;
A2=find(tA>0 & tA<tmax+1);
tAm((tA(A2)-1)*n+A2)=1;
tAm=sparse(tAm);
% figure; spy(tEm); hold on; spy(tAm,'r'); axis square
if inclLST
    % LST conversion
    tLm=false(n,tmax);
    tLm((tL(L)-1)*n+L)=1;
end
% % Relapse
% tRLm=false(n,tmax);
% tRLm((tRL(RL)-1)*n+RL)=1;
% % Relapse recovery
% tRLRm=false(n,tmax);
% tRLRm((tRLR(RLR)-1)*n+RLR)=1;
% PKDL
tPm=false(n,tmax);
tPm((tP(P)-1)*n+P)=1;
% Birth
tBm=false(n,tmax);
tBm((tB(B)-1)*n+B)=1;
% Death
tDm=false(n,tmax);
tDm((tD(D)-1)*n+D)=1;
% Pre-birth
preB=false(n,tmax);
preB(B,:)=1-cumsum(tBm(B,:),2);
% Immigration
tIMm=false(n,tmax);
tIMm((tIM(IM)-1)*n+IM)=1;
% Pre-immigration
preIM=false(n,tmax);
preIM(IM,:)=1-cumsum(tIMm(IM,:),2);
% Emigration
tEMm=false(n,tmax);
tEMm((tEM(EM)-1)*n+EM)=1;
% % Entry
% tENm=false(n,tmax);
% tENm((tEN(EN)-1)*n+EN)=1;
% % Exit
% tEXm=false(n,tmax);
% tEXm((tEX(EX)-1)*n+EX)=1;

% Make susceptibility status matrix
if inclLST
%     S=cumsum(tENm,2)-max(max(max(cumsum(tEm,2),cumsum(tPm,2)),cumsum(tLm,2)),cumsum(tEXm,2)); % remove LST+ individuals from susceptibles
    S=1-max(preB,preIM)-max(max(max(max(max(cumsum(tEm,2),cumsum(tAm,2)),cumsum(tPm,2)),cumsum(tDm,2)),cumsum(tLm,2)),cumsum(tEMm,2)); % remove LST+ individuals from susceptibles
else
%     S=cumsum(tENm,2)-max(max(cumsum(tEm,2),cumsum(tPm,2)),cumsum(tEXm,2)); % don't remove LST+ individuals
    S=1-max(preB,preIM)-max(max(max(max(cumsum(tEm,2),cumsum(tAm,2)),cumsum(tPm,2)),cumsum(tDm,2)),cumsum(tEMm,2)); % don't remove LST+ individuals
end
S(prevK,:)=0; % remove previous KA cases from susceptibles
S(prevA,:)=0; % remove previous asymptomatic infections from susceptibles
S(actvA,:)=0; % remove initially active asymptomatic infections from susceptibles
S(tI<=maxIP,:)=0; % remove cases with onset before maxIP (also excludes cases with active KA at start of study and a couple of cases with onset in 2002 before immigration)
S(IpreEXTIM,:)=0; % remove KA cases with onset before or at migration in
S(EXTIMsoonI,:)=0; % remove KA cases with onset within 6 months of migration in
S(IpreINTIM,:)=0; % remove KA cases with onset before internal migration in
S(PpreINTIM,:)=0; % remove PKDL cases with onset before internal migration in
S(PpreEXTIM,:)=0; % remove PKDL cases (without prior KA) with onset before external migration in
S(prevAPA,:)=0; % remove previous asymptomatic infections from susceptibles
S(actvAPA,:)=0; % remove initially active asymptomatic infections from susceptibles
S(IM_INprevAactvA,:)=0; % remove susceptibility from 2nd observations of internal migrators asymptomatically infected before the start of the study
S(IM_IN(ismember(IM_OUT,Asx)),:)=0; % remove susceptibility from 2nd observations of internal migrators asymptomatically infected during 1st observation

% figure; histogram(S)

% % Checks
% figure; spy(S(setdiff(I,IandP),:)); axis square; pause(2); hold on
% spy(h(~ismember(I,IandP),:),'r'); axis square; hold off
% figure; spy(S(IandP,:)); axis square; pause(2); hold on
% spy(h(ismember(I,IandP),:),'r'); axis square; hold off
% figure; spy(S([IMI-1,IMI,IMP(1)-1,IMP(1),IMP(2)-1,IMP(2)],:)); axis square; pause(2); hold on
% spy(h([find(I==IMI-1),end-2,find(I==IMP(1)-1),end-1,nI+find(PNI==IMP(2)-1),end],:),'r'); axis square; hold off

% Susceptibility vector
% Variation in susceptibility according to age
% s=load('~/Dropbox/Visceral Leishmaniasis/CarynBernData/2004data/SpatiotemporalModelling/data_final.mat');
% [pars2,~]=FitCatModLST2(s.data);
% age=max((1-tB)/12,0);
% sigma=1-pars2(2)*(1-exp(-pars2(1)*age));
% S=bsxfun(@times,S,sigma);
% % No variation in susceptibility
% sigma=ones(n,1);

S0=intersect(union(Sus,find(S(:,1))),Pres0);

%% SET UP MCMC
burnin=round(niters/10); %1e3; %1e5; %1e4; % burn-in
nu=numel(u); % number of parameters to update
% pname={'beta','alpha','epsilon','delta'}; % parameter labels
% pname={'beta','alpha','epsilon','delta','lambda_0'};
pname={'beta','alpha','epsilon','delta','lambda_0','h_4','p_I'};
np=numel(pname); % number of possible parameters in model

% PRIORS
% CHANGE PRIOR MEAN FOR alpha TO 100m !!
% prior_mean=[1,50,1,1]; % prior distribution means
prior_shape=ones(1,6);%[1,1,1,1,10,1];%[1,1,1,1,1];%[1,1,1,1];%[1,2,1,1];
prior_scale=[10,100,1,1,0.01,h40];%[10,100,1,1,lambda0/prior_shape(5),h40];%[1,100,1,1,lambda0/prior_shape(5),h40];%[1,100,1,1,lambda0/prior_shape(5),1];%[1,100,1,1,pars/(12*prior_shape(5))];%[1,100,1,1,1];%[1,100,1,1];%[1,50,1e-3,0.1];
priorpdf={'gam','gam','gam','gam','gam','gam','beta'};
priorp=cell(1,np);
for i=1:6
    priorp{i}=[prior_shape(i),prior_scale(i)];
end
% priorp{7}=[15,85];
priorp{7}=[1,1];

% PROPOSALS
% ppmean=zeros(1,np); % vector for storing means for proposal distributions for transmission parameters
% ppmean=[beta0,alpha0,epsilon0,delta0,lambda0];%[beta0,alpha0,epsilon0,delta0];
% ppvar0=diag([beta0,alpha0,epsilon0,delta0,lambda0].^2);%diag([beta0,alpha0,epsilon0,delta0].^2); % initial variances for proposal distribution for block update
% ppvar0=diag([0.03,1/3,1,1,1].*[beta0,alpha0,epsilon0,delta0,lambda0].^2);%diag([beta0,alpha0,epsilon0,delta0].^2); % initial variances for proposal distribution for block update
% ppvar0=diag([0.03,1/3,1,1,1].*[0.2,100,1e-4,0,lambda0].^2);%diag([beta0,alpha0,epsilon0,delta0].^2); % initial variances for proposal distribution for block update
% ppvar0=diag([0.2,100,1e-4,0,lambda0].^2);

% ppvar0=diag([0.2,2e3,1e-6,0,lambda0^2]);
% ppvar0=diag([0.01,2e3,1e-6,0,lambda0^2]);
% ppvar0=diag([0.01,400,1e-7,0,0]);
% ppvar0=diag([0.01,100,1e-8,0,0,1e-4]);
% ppvar0=diag([0.01,100,1e-8,0,0,1e-4,1e-4]);
% ppvar0=diag([0.01,100,1e-7,0,0,1e-4,1e-4]);
% ppvar0=diag([0.01,100,1e-7,1e-4,0,1e-4,1e-4]);
% % Proposal variances for para 1
% ppvar0=diag([0.01,100,1e-7,4e-4,0,1e-4,1e-4]);
% Proposal variances for full data
ppvar0=diag([0.01,100,1e-7,4e-4,4e-8,1e-4,1e-4]/25);
% ppvar0=diag([0.01,100,1e-8,0,0,1e-4,1e-3]);
% ppvar0=0.1^2*diag([beta0,alpha0,epsilon0,delta0,lambda0].^2)/nu;

% ppvar1=diag([beta0,alpha0,epsilon0,delta0].^2); % initial variances for proposeal distributions for single component updates
% ppvar=zeros(np); % vector for storing variances for proposal distributions for transmission parameters
ppvar=ppvar0; % vector for storing variances for proposal distributions for transmission parameters
c=NaN(niters+1,1); % initial scale factor for proposal covariance matrix for block update 
c(1)=1;
ppvar1=0.2;
c1=NaN(niters+1,1);
c1(1)=1;
nEmoves=round(nOR/5); % number of cases with onset and recovery times to propose new infection times for in each step 
pick=zeros(nEmoves+nNONR+nONR+nRNO,niters);
for i=1:niters 
pick(:,i)=[OR(randperm(nOR,nEmoves));NONR;ONR;RNO]; % indices of individuals to propose new infection times for in each step
% N.B. randperm rather than randi to avoid possibility of changing same infection time twice in one step
end
% pick=[OR(randi(nOR,nEmoves,niters));repmat(NONR,1,niters);repmat(ONR,1,niters);repmat(RNO,1,niters)];
ERvar=4; % variance for infection period moves
nAupdts=round(nAsx/5);%round(nAsx/10);%round((round((1-pI)/pI*nI)-nPA)/5);%round(nAsx/3);%round(nAsx/10);%round(nAsx/20);%
pickA=zeros(nAupdts,niters,'uint16');
% for i=1:niters
% pickA(:,i)=SusA(randperm(nSusA,nAupdts));
% end
% % pickA=SusA(randi(nSusA,nAupdts,niters));
M=tmax; %2; %0; % maximum allowable asymptomatic infection time jump

% eta=1; %0.95; % probability of random uniform vs infectious pressure proposal for asymptomatic infection

% INITIALISATION
% Fitting 4 (or 5) transmission parameters: beta, alpha, epsilon (and delta) and p1
% pold=[beta0,alpha0,epsilon0,delta0]; % set old parameter values
% pold=[beta0,alpha0,epsilon0,delta0,lambda0];
% pold=[beta0,alpha0,epsilon0,delta0,0.001/12];
% pold=[beta0,alpha0,epsilon0,delta0,lambda0,h40];
pold=[beta0,alpha0,epsilon0,delta0,lambda0,h40,pI0];
ppmean=pold;
% plb=zeros(1,np); % lower limits for parameters
% pub=Inf(1,np); % upper limits for parameters
plb=[zeros(1,6),0]; % lower limits for parameters
pub=[Inf(1,6),1]; % upper limits for parameters
p1new=p10; % initial p1 value
LL1old=L1(S,lambda);
% LL2old=L2(lambda,pI,tEm);
LL2old=L2(lambda,pI0,tEm);
LL3old=L3(IPold(I),r1,p10);
% LL4old=L4(lambda,pI,tAm);
LL4old=L4(lambda,pI0,tAm);
LL5old=L5(age,S0,actvA,prevA,lambda0,pI0,p2);
LL6old=L5(age,S0PA,actvAPA,prevAPA,lambda0,pI0,p2);
 
%%
% Matrices and vectors for saving parameter and log-likelihood values
% p=zeros(niters,np); % matrix for saving parameter values
p=zeros(niters+1,np); % matrix for saving parameter values
p(1,:)=pold;
p1=zeros(niters,1); % vector for saving p1 values
K0=zeros(niters,1); % vector for saving spatial kernel normalisation constants
LL=zeros(niters,1); % vector for saving log-likelihood values
% terms=zeros(niters,5); % matrix for saving individual log-likelihood terms
terms=zeros(niters,6); % matrix for saving individual log-likelihood terms
% tEs=zeros(nI,niters);
tEs=zeros(nI,niters,'int8'); % integer matrix to save memory
% tAs=NaN(n,niters);
% tRAs=NaN(n,niters);
% Matrices for saving asymptomatic infection and recovery times ? need to 
% be integer matrices as otherwise use too much memory
tAs=zeros(n,niters,'uint8'); 
tAs(setdiff((1:n)',[SusA;PA]),:)=tmax+2; % use dummy asymptomatic infection time for symptomatic individuals
tRAs=zeros(n,niters,'uint8');
tRAs(setdiff((1:n)',[SusA;PA]),:)=tmax+2; % use dummy asymptomatic recovery time for symptomatic individuals
tIsNONR=zeros(nNONR,niters);
tRsNONR=zeros(nNONR,niters);
tIsRNO=zeros(nRNO,niters);
tRsONR=zeros(nONR,niters);
tIsANONR=zeros(nANONR,niters);
tRsANONR=zeros(nANONR,niters);
tRsAONR=zeros(nAONR,niters);
tRLsRLO=zeros(nRLO,niters); % this is unnecessary as these values don't change!
tRLRsRLO=zeros(nRLO,niters);
tRLsRLNO=zeros(nRLNO,niters);
tRLRsRLNO=zeros(nRLNO,niters);
 
% Initialise new status and infectiousness matrices
pnew=pold;
rateHH_new=rateHH;

rateHHA_new=rateHHA;
rateHHPA_new=rateHHPA;

lambdaHH_new=lambdaHH;
lambdaHHA_new=lambdaHHA;
lambdaHHPA_new=lambdaHHPA;
lambda_new=lambda;

% lambdaHH_new1=lambdaHH_new;

lambda_mean=lambda;

tEnew=tE;
tAnew=tA;
tRAnew=tRA;
tInew=tI;
tRnew=tR;
tRLnew=tRL;
tRLRnew=tRLR;
tEmnew=tEm;
tAmnew=tAm;
Snew=S; 
IPnew=IPold;
hnew=h;
hAnew=hA;
hPAnew=hPA;

% hA1new=hA1;

A1new=A1;
nA1new=nA1;
actvAnew=actvA;
prevAnew=prevA;
Susnew=Sus;
S0new=S0;
actvAPAnew=actvAPA;
prevAPAnew=prevAPA;
S0PAnew=S0PA;
RAobs2actvAnew=RAobs2actvA;
RAobs2new=RAobs2;

% Initialise acceptance and rejection counts and rates
acc_p=0; % initial acceptance count for transmission parameter updates
rej_p=0; % initial rejection count for transmission parameter updates
acc_rate_p=0; % initial acceptance rate for transmission parameter updates
acc_E=0;
rej_E=0;
acc_rate_E=0;
acc_Aadd=0;
rej_Aadd=0;
acc_rate_Aadd=0;
acc_Arem=0;
rej_Arem=0;
acc_rate_Arem=0;
acc_Amov=0;
rej_Amov=0;
acc_add=zeros(n,1);
rej_add=zeros(n,1);
acc_rem=zeros(n,1);
rej_rem=zeros(n,1);
acc_mov=zeros(n,1);
rej_mov=zeros(n,1);
acc_rate_Amov=0;
acc_I=0;
rej_I=0;
acc_rate_I=0;
acc_ERmove=0;
rej_ERmove=0;
acc_rate_ERmove=0;
acc_R=0;
rej_R=0;
acc_rate_R=0;
acc_AIRmove=0;
rej_AIRmove=0;
acc_rate_AIRmove=0;
acc_AR=0;
rej_AR=0;
acc_rate_AR=0;
acc_PA=0;
rej_PA=0;
acc_rate_PA=0;
acc_RLNO=0;
rej_RLNO=0;
acc_rate_RLNO=0;
acc_RLO=0;
rej_RLO=0;
acc_rate_RLO=0;
% acc_beta=0;
% rej_beta=0;
% acc_rate_beta=0;

nbins=50;
scrnsz=get(0,'ScreenSize');

% error_lambda=NaN(niters,1);

% fwd=0:tmax+1;

% % Get coefficients for linear relationship between mean asymptomatic 
% % infection time and beta
% d=FitLinMod('MCMC_NBIP_PKDL_ASX110.mat');

%% MCMC LOOP
for k=1:niters
%     k
    %% UPDATE TRANSMISSION PARAMETERS USING ADAPTIVE RANDOM WALK METROPOLIS-HASTINGS    
%     if k<=max(burnin,100) || rand<0.05 %k<=2*burnin || rand<0.05
% %         pnew(u)=mvnrnd(pold(u),0.1^2*ppvar0(u,u)/nu);
%         pnew(u)=mvnrnd(pold(u),ppvar0(u,u));
%     else
        pnew(u)=mvnrnd(pold(u),c(k)^2*2.38^2*ppvar(u,u)/nu);
%         pnew(u)=mvnrnd(pold(u),c(k)^2*2.38^2*diag(ppvar(u,u))'/nu);
%     end
    
    if all(pnew(u)>plb(u) & pnew(u)<pub(u)) % check if prior probability is non-zero before calculating log-likelihood
        % Calculate ratio of prior probabilities for new and old value
%         q=sum(log(gampdf(pnew(u),prior_shape(u),prior_scale(u)))-log(gampdf(pold(u),prior_shape(u),prior_scale(u))));
        q=0;
        for i=u
            if numel(priorp{i})==1
                q=q+log(pdf(priorpdf{i},pnew(i),priorp{i}))-log(pdf(priorpdf{i},pold(i),priorp{i}));
            elseif numel(priorp{i})==2
                q=q+log(pdf(priorpdf{i},pnew(i),priorp{i}(1),priorp{i}(2)))-log(pdf(priorpdf{i},pold(i),priorp{i}(1),priorp{i}(2)));
            end
        end
        
        [KHH_new,K0new]=Knl_fast(dHH,dHHsqrd,pnew(2),pnew(1),typ,n,nHH,ib,f,1:n);
        rateHHA_new=pnew(1)*KHH_new;
        if ismember(4,u)
            rateHHA_new=rateHHA_new+pnew(4)*d0;
        end
        rateHH_new=rateHHA_new(:,ib(IPNIA));
        rateHHPA_new=rateHHA_new(:,ib(PA));
        if ismember(6,u)
            hAnew=pnew(6)/pold(6)*hA;
        end
        lambdaHHA_new=rateHHA_new(:,ib([A1;RAobs2actvA;RAobs2]))*hAnew([A1;RAobs2actvA;RAobs2],:);
        lambdaHHPA_new=rateHHPA_new*hPA;
        lambdaHH_new=rateHH_new*h+lambdaHHPA_new+lambdaHHA_new+pnew(3); % update whole infectious pressure
        lambda_new=lambdaHH_new(ib,:); % expand infectious pressure

        LL1new=L1(S,lambda_new);
%         LL2new=L2(lambda_new,pI,tEm);
%         LL4new=L4(lambda_new,pI,tAm);
        LL2new=L2(lambda_new,pnew(7),tEm);
        LL4new=L4(lambda_new,pnew(7),tAm);
        LL5new=L5(age,S0,actvA,prevA,pnew(5),pnew(7),p2);
        LL6new=L5(age,S0PA,actvAPA,prevAPA,pnew(5),pnew(7),p2);
        LLnew=LL1new+LL2new+LL4new+LL5new+LL6new;
        LLold=LL1old+LL2old+LL4old+LL5old+LL6old;

        log_ap=LLnew-LLold+q; % calculate Metropolis-Hastings acceptance probability

        if log_ap > log(rand) % for acc_prob<=1, change if greater than rand
            pold=pnew; % keep updated parameter values
            lambda=lambda_new; lambdaHH=lambdaHH_new; rateHH=rateHH_new; rateHHA=rateHHA_new; rateHHPA=rateHHPA_new; KHH=KHH_new; K0old=K0new;
            lambdaHHA=lambdaHHA_new; lambdaHHPA=lambdaHHPA_new;
            hA=hAnew;
            LL1old=LL1new; LL2old=LL2new; %LL3old=LL3new; 
            LL4old=LL4new; LL5old=LL5new; LL6old=LL6new; % keep updated info
            acc_p=acc_p+1; % add to acceptance
%             c(k+1)=1;
%             c(k+1)=(1+200/(1000+k))*c(k);
            c(k+1)=(1+100/(100+k))*c(k);
        else
            lambda_new=lambda; lambdaHH_new=lambdaHH; %rateHH_new=rateHH; rateHHA_new=rateHHA; % revert lambda_new to old value as only part of lambda_new is updated in infection time updates below
            lambdaHHA_new=lambdaHHA; lambdaHHPA_new=lambdaHHPA;
            hAnew=hA;
            rej_p=rej_p+1; % add to rejection
%             c(k+1)=1;
%             c(k+1)=(1+200/(1000+k))^(0.234/(0.234-1))*c(k);
            c(k+1)=(1+100/(100+k))^(0.234/(0.234-1))*c(k);
%             c(k+1)=max(1,(1+200/(1000+k))^(0.234/(0.234-1))*c(k));
%             c(k+1)=max(1,(1+100/(100+k))^(0.234/(0.234-1))*c(k));
        end
    else % if prior probability is 0, reject immediately
        rej_p=rej_p+1;
%         c(k+1)=1;
%         c(k+1)=(1+200/(1000+k))^(0.234/(0.234-1))*c(k);
        c(k+1)=(1+100/(100+k))^(0.234/(0.234-1))*c(k);
%         c(k+1)=max(1,(1+200/(1000+k))^(0.234/(0.234-1))*c(k));
%         c(k+1)=max(1,(1+100/(100+k))^(0.234/(0.234-1))*c(k));%
    end
    
    %% GIBBS UPDATE SUCCESS PROBABILITY PARAMETER FOR ASYMPTOMATIC INFECTION PERIOD DISTRIBUTION
    p1new=betarnd(a+r1*nI,b+sum(IPold(I))-nI);
 
    %% UPDATE INFECTION TIMES    
    for i=1:size(pick,1)
        j=pick(i,k); % get index of infection time to update
        tEnew(j)=tI(j)-(nbinrnd(r1,p1new)+1); % draw new incubation period from negative binomial distribution
        
        if tEnew(j)>=tB(j)+1 && ~(tI(j)>maxIP && tEnew(j)<1) && ~(tI(j)-tIM(j)-1>maxIP_IM && tEnew(j)<tIM(j)+1) % calculate log-likelihood if new infection time is after birth, not before start of study if onset is after max_IP, and not before immigration if onset is more than maxIP_IM months after immigration
            tEj=tE(j);
            tEjnew=tEnew(j);
            tIMj=tIM(j);
            IPnew(j)=tI(j)-tEjnew; % new incubation period
            q=log(nbinpdf(IPold(j)-1,r1,p1new))-log(nbinpdf(IPnew(j)-1,r1,p1new)); % proposal ratio for infection time update
            
            if tI(j)>maxIP && ~IpreEXTIM(j) && ~EXTIMsoonI(j) % only update infection time and susceptibility matrices for cases with onset after initial window and after entry
                tEmnew(j,tEj)=0; % remove old infection time
                tEmnew(j,tEjnew)=1; % add new infection time
                Snew(j,tEjnew:tEj-1)=0; % remove old susceptible times if new infection time is earlier (from month indvdl is infctd up to but not incl. month before old infection time)
%                 Snew(j,tEj:tEjnew-1)=sigma(j); % add new susceptible times if new infection time is later (up to but not incl. month indvdl is infctd)
                Snew(j,tEj:tEjnew-1)=1; % add new susceptible times if new infection time is later (up to but not incl. month indvdl is infctd)
            end
            
            m=(IPNIA==j);%ismember(IPNIA,j);
            entry=max(0,tIMj);
            hnew(m,max(entry,tEj)+1:tEjnew)=0; % remove old infectiousness if new infection time is later
            hnew(m,max(entry,tEjnew)+1:tEj)=h0; % add infectiousness if new infection time is earlier
            erlrE=max(entry,min(tEj,tEjnew)); % index of column for earlier infctn time between old and new infection time
            ltrE=max(entry,max(tEj,tEjnew)); % index of column for later infctn time between old and new infection time
            idx=erlrE+1:ltrE;
            lambdaHH_new(:,idx)=rateHH*hnew(:,idx)+lambdaHHPA(:,idx)+lambdaHHA(:,idx)+pold(3); % update infectious pressure
            lambda_new(:,idx)=lambdaHH_new(ib,idx); % expand infectious pressure

            LL1new=L1(Snew,lambda_new);
%             LL2new=L2(lambda_new,pI,tEmnew);
            LL2new=L2(lambda_new,pold(7),tEmnew);
            LL3new=L3(IPnew(I),r1,p1new);
%             LL4new=L4(lambda_new,pI,tAm);
            LL4new=L4(lambda_new,pold(7),tAm);
            LLnew=LL1new+LL2new+LL3new+LL4new;
            LLold=LL1old+LL2old+LL3old+LL4old;
            
            log_ap=LLnew-LLold+q; % calculate Metropolis-Hastings acceptance probability
            
            if log_ap > log(rand) % for acc_prob<=1, change if greater than rand
                IPold(j)=IPnew(j); tE(j)=tEnew(j); S(j,:)=Snew(j,:); tEm(j,:)=tEmnew(j,:);
                h(m,:)=hnew(m,:); lambda(:,idx)=lambda_new(:,idx); lambdaHH(:,idx)=lambdaHH_new(:,idx);
                LL1old=LL1new; LL2old=LL2new; LL3old=LL3new; 
                LL4old=LL4new; % keep updated info
                acc_E=acc_E+1;
            else
                IPnew(j)=IPold(j); tEnew(j)=tE(j); Snew(j,:)=S(j,:); tEmnew(j,:)=tEm(j,:);
                hnew(m,:)=h(m,:); lambda_new(:,idx)=lambda(:,idx); lambdaHH_new(:,idx)=lambdaHH(:,idx); % keep old values, don't change log-likelihood
                rej_E=rej_E+1;
            end
        else % otherwise reject immediately
            tEnew(j)=tE(j);
            rej_E=rej_E+1;
        end
    end
     
    %% UPDATE ASYMPTOMATIC INFECTION TIMES    
    % Recalculate cumulative infectious pressure
%     Lambda=sum(lambda_mean.*rngm,2);
    CumLambda=cumsum(lambda_mean.*rngm,2);
    CumLambda=[zeros(n,1),CumLambda(:,1:end-1)];
    probA1=[exp(-CumLambda).*(1-exp(-(1-pold(7))*lambda_mean.*rngm)),exp(-sum(lambda_mean.*rngm,2))];
    probA1=bsxfun(@rdivide,probA1,sum(probA1,2));
    prob0=ProbInitStatus(age,pold(5),p2);
    prob0=[prob0,1-sum(prob0,2)];
    probA=[1-prob0(:,1),bsxfun(@times,prob0(:,1),probA1)];
    Lambda=sum(probA(:,1:end-1),2);
    for i=1:nAupdts
% %        j=SusA(randi(nSusA));
%        j=pickA(i,k);
% %        j=SusA(randsample(nSusA,1,true,LambdaI(SusA)));
%        j=SusA(randsample(nSusA,1,true,Lambda(SusA)));
       j=SusA(sum(rand>=cumsum(Lambda(SusA))/sum(Lambda(SusA)))+1);
       mig_out=ismember(j,IM_OUT);
       mig_in=ismember(j,IM_IN);
       pickA(i,k)=j;
       if ~(mig_out || mig_in) || (mig_out && ~ismember(IM_IN(IM_OUT==j),[prevA;A1])) || (mig_in && ~ismember(IM_OUT(IM_IN==j),[prevA;A1]))
       t=tA(j);
       s=tRA(j);
       idx1=t+1:min(min(s,rng(j,2)),tmax);
       if mig_out 
           j1=IM_IN(IM_OUT==j);
           if t<=rng(j,2)-1 && s>rng(j,2)-1 % check if individual was infected before internal migration and remained infected after moving
               idx3=rng(j,2)+1:min(min(s,rng(j1,2)),tmax); % create index for infection during 2nd observation
           end
%        elseif mig_in
%            j2=IM_OUT(IM_IN==j);
       end
% %        tp=randsample(0:tmax+1,1,true,probA(j,:)); %t; %
%        tp=fwd(sum(rand>=cprobA(j,:))+1);
       if t==0 % currently asymptomatically infected before start of study (N.B. means must be individual's 1st observation)
%            fwd=[1:rng(j,2)-1,tmax+1];
% %            tp=fwd(randi(numel(fwd)));
%            tp=fwd(randsample(numel(fwd),1,true,probA(j,fwd+1)));
%            fwd=0:tmax+1;
%            tp=fwd(sum(rand>=cprobA(j,:))+1);
           fwd=[0:rng(j,2)-1,tmax+1];
           tp=fwd(sum(rand>=cumsum(probA(j,fwd+1))/sum(probA(j,fwd+1)))+1);
           if tp==0
%                bck=0:tmax+1;
               bck=[0:rng(j,2)-1,tmax+1];
               if rand<prob0(j,3)/sum(prob0(j,2:3)) % recovered before t=0 
                   tRAnew(j)=0;
                   if s==0 % if individual is not already in asymptomatic infectiousness matrix as they recovered by t=0
                       idx=[];
                       q=0;
                   else % if individual is already in asymptomatic infectiousness matrix as they recovered after t=0
                       prevAnew=[prevA;j];
                       actvAnew=actvA(actvA~=j);
                       A1new=A1(A1~=j);
                       nA1new=nA1-1;
                       idx=idx1;
                       hAnew(j,idx1)=0;
                       % CHECK INDICES!!                       
                       if mig_out && s>rng(j,2)-1 % still need to check if individual migrated!
                           RAobs2actvAnew=RAobs2actvA(RAobs2actvA~=j1); % remove from set of 2nd observations with infectiousness
                           hAnew(j1,idx3)=0; % remove infectiousness if currently still infected after moving
                           idx=[idx,idx3];
                           tmp=rateHHA(:,ib(j))*hA(j,idx)+rateHHA(:,ib(j1))*hA(j1,idx);
                       else
                           tmp=rateHHA(:,ib(j))*hA(j,idx);
                       end
                       lambdaHH_new(:,idx)=lambdaHH(:,idx)-tmp;
                       lambdaHHA_new(:,idx)=lambdaHHA(:,idx)-tmp;
                       lambda_new(:,idx)=lambdaHH_new(ib,idx);
                       q=log(prob0(j,2)/sum(prob0(j,2:3)))-log(prob0(j,3)/sum(prob0(j,2:3)));
                   end
               else % recovered after t=0
                   if mig_out
                       probAIPnew=[geopdf(0:rng(j1,2)-2,p2),1-geocdf(rng(j1,2)-2,p2)];
                       fwdRA=[1:rng(j1,2)-1,tmax+1];                       
                   else
                       probAIPnew=[geopdf(0:rng(j,2)-2,p2),1-geocdf(rng(j,2)-2,p2)];                       
                       fwdRA=[1:rng(j,2)-1,tmax+1];
                   end
%                    tRAnew(j)=fwdRA(randsample(numel(fwdRA),1,true,probAIPnew));
                   tRAnew(j)=fwdRA(sum(rand>=cumsum(probAIPnew))+1);
%                    tRAnew(j)=durA;
                   if s==0 % if individual is not already in asymptomatic infectiousness matrix as they recovered by t=0
%                       if ismember(j,IM_OUT) && tRAnew(j)>rng(j,1)-1
%                           idx=1:min(min(tRAnew(j),rng(j1,2)),tmax);
%                           hAnew(j,1:rng(j,2)-1)=pold(6);
%                           hAnew(j1,rng(j,2):idx(end))=pold(6);
%                           tmp=rateHHA(:,j)*hAnew(j,idx);
%                       else
%                           idx=1:min(min(tRAnew(j),rng(j,2)),tmax);
%                           hAnew(j,idx)=pold(6);
%                           tmp=rateHHA(:,j)*hAnew(j,idx);
%                       end
                      prevAnew=prevA(prevA~=j);
                      actvAnew=[actvA;j];
                      A1new=[A1;j];
                      nA1new=nA1+1;
                      idx=1:min(min(tRAnew(j),rng(j,2)),tmax);
% %                       hAnew(j,idx)=h4;
                      hAnew(j,idx)=pold(6);
                      if mig_out && tRAnew(j)>rng(j,2)-1
                          RAobs2actvAnew=[RAobs2actvA;j1];
                          idx4=rng(j,2)+1:min(min(tRAnew(j),rng(j1,2)),tmax);
                          hAnew(j1,idx4)=pold(6);
                          idx=[idx,idx4];
                          tmp=rateHHA(:,ib(j))*hAnew(j,idx)+rateHHA(:,ib(j1))*hAnew(j1,idx);
                      else
                          tmp=rateHHA(:,ib(j))*hAnew(j,idx);
                      end
                      lambdaHH_new(:,idx)=lambdaHH(:,idx)+tmp;
                      lambdaHHA_new(:,idx)=lambdaHHA(:,idx)+tmp;
                      lambda_new(:,idx)=lambdaHH_new(ib,idx);
                      q=log(prob0(j,3)/sum(prob0(j,2:3)))-log(prob0(j,2)/sum(prob0(j,2:3)));%+log(probAIPnew(tRAnew(j)));                 
                  else % individual is already in asymptomatic infectiousness matrix as they recovered after t=0
                      idx2=1:min(min(tRAnew(j),rng(j,2)),tmax);
                      hAnew(j,idx1)=0;
%                       hAnew(j,idx2)=h4;
                      hAnew(j,idx2)=pold(6);
                      idx=union(idx1,idx2);
                      
                      if ~mig_out || (mig_out && s<=rng(j,2)-1 && tRAnew(j)<=rng(j,2)-1)
                          tmp=rateHHA(:,ib(j))*(hAnew(j,idx)-hA(j,idx));
                      else % individual internally migrated and old or new recovery time is after migration time
                          if s<=rng(j,2)-1 && tRAnew(j)>rng(j,2)-1
                              RAobs2actvAnew=[RAobs2actvA;j1];
                              idx4=rng(j,2)+1:min(min(tRAnew(j),rng(j1,2)),tmax);
                              hAnew(j1,idx4)=pold(6);
                              idx=[idx,idx4];
                              tmp=rateHHA(:,ib(j))*(hAnew(j,idx)-hA(j,idx))+rateHHA(:,ib(j1))*hAnew(j1,idx);
                          elseif s>rng(j,2)-1 && tRAnew(j)<=rng(j,2)-1
                              RAobs2actvAnew=RAobs2actvA(RAobs2actvA~=j1);
                              hAnew(j1,idx3)=0; % remove infectiousness if currently still infected after moving
                              idx=[idx,idx3];
                              tmp=rateHHA(:,ib(j))*(hAnew(j,idx)-hA(j,idx))-rateHHA(:,ib(j1))*hA(j1,idx);
                          else %if s>rng(j,2)-1 && tRAnew(j)>rng(j,2)-1
                              hAnew(j1,idx3)=0; % remove infectiousness if currently still infected after moving
                              idx4=rng(j,2)+1:min(min(tRAnew(j),rng(j1,2)),tmax);
                              hAnew(j1,idx4)=pold(6);
                              idx=[idx,union(idx3,idx4)];
                              tmp=rateHHA(:,ib(j))*(hAnew(j,idx)-hA(j,idx))+rateHHA(:,ib(j1))*(hAnew(j1,idx)-hA(j1,idx));
                          end
                      end
                      
                      lambdaHH_new(:,idx)=lambdaHH(:,idx)+tmp;
                      lambdaHHA_new(:,idx)=lambdaHHA(:,idx)+tmp;
                      lambda_new(:,idx)=lambdaHH_new(ib,idx);
                      q=0; % backward and forward proposal probabilities cancel
                  end              
               end
           elseif tp==tmax+1 % proposed not asymptomatically infected before or during study (tp=tmax+1)
               tRAnew(j)=tmax+1;
%                bck=0:rng(j,2)-1;
%                bck=0:tmax+1;
               bck=[0:rng(j,2)-1,tmax+1];
               if s==0 % if individual is not already in asymptomatic infectiousness matrix as they recovered by t=0
                   prevAnew=prevA(prevA~=j);
                   idx=[];
                   q=log(prob0(j,3)/sum(prob0(j,2:3)));
               else % individual is already in asymptomatic infectiousness matrix as they recovered after t=0
%                    probAIP=[geopdf(0:tmax-1,p2),1-geocdf(tmax-1,p2)];
                   actvAnew=actvA(actvA~=j); % remove from initially actively asymptomatic individuals
                   A1new=A1(A1~=j);
                   nA1new=nA1-1;
                   idx=idx1;
                   hAnew(j,idx)=0; % remove row for previous infectiousness
                   if mig_out && s>rng(j,2)-1 % still need to check if individual migrated!
                       RAobs2actvAnew=RAobs2actvA(RAobs2actvA~=j1);
                       hAnew(j1,idx3)=0; % remove infectiousness if currently still infected after moving
                       idx=[idx,idx3];
                       tmp=rateHHA(:,ib(j))*hA(j,idx)+rateHHA(:,ib(j1))*hA(j1,idx);
                   else
                       tmp=rateHHA(:,ib(j))*hA(j,idx);
                   end
                   lambdaHH_new(:,idx)=lambdaHH(:,idx)-tmp;
                   lambdaHHA_new(:,idx)=lambdaHHA(:,idx)-tmp; 
                   lambda_new(:,idx)=lambdaHH_new(ib,idx);
                   q=log(prob0(j,2)/sum(prob0(j,2:3))); %+log(probAIP(s));
               end
               Susnew=[Sus;j]; %sort([Sus;j]); % add j to intially susceptible set
               S0new=[S0;j];
               if mig_out
                   Snew(j1,rng(j,2):rng(j1,2)-1)=1; % add susceptibility to 2nd observation now that individual is not asymptomatically infected before start of study
               end
           else % proposed asymptomatically infected during study (tp in [1,tmax])
               if mig_out
                   probAIPnew=[geopdf(0:rng(j1,2)-tp-2,p2),1-geocdf(rng(j1,2)-tp-2,p2)];
                   fwdRA=[tp+1:rng(j1,2)-1,tmax+1];
               else
                   probAIPnew=[geopdf(0:rng(j,2)-tp-2,p2),1-geocdf(rng(j,2)-tp-2,p2)];
                   fwdRA=[tp+1:rng(j,2)-1,tmax+1];
               end
%                tRAnew(j)=fwdRA(randsample(numel(fwdRA),1,true,probAIPnew));
%                tRAnew(j)=min(tp+durA,tmax+1);
               tRAnew(j)=fwdRA(sum(rand>=cumsum(probAIPnew))+1);
               idx2=tp+1:min(min(tRAnew(j),rng(j,2)),tmax);
               if mig_out && tRAnew(j)>rng(j,2)-1
                   idx4=rng(j,2)+1:min(min(tRAnew(j),rng(j1,2)),tmax);
               end
%                bck=[0,max(1,tp-M):tp-1,tp+1:min(rng(j,2),tp+M),tmax+1];
%                bck=[0:tp-1,tp+1:rng(j,2)-1,tmax+1];
               bck=[0,max(1,tp-M):min(rng(j,2)-1,tp+M),tmax+1];
               tAmnew(j,tp)=1; % add new infection time
               if s==0 % if individual is not already in asymptomatic infectiousness matrix as they recovered by t=0
                   prevAnew=prevA(prevA~=j); % remove j from previously asymptomatically infected set
                   A1new=[A1;j];
                   nA1new=nA1+1;
                   idx=idx2;
%                    hAnew(j,idx)=h4; % add new row for their infectiousness
                   hAnew(j,idx)=pold(6); % add new row for their infectiousness
                   if mig_out && tRAnew(j)>rng(j,2)-1
                       RAobs2new=[RAobs2;j1];
                       hAnew(j1,idx4)=pold(6);
                       idx=[idx,idx4];
                       tmp=rateHHA(:,ib(j))*hAnew(j,idx)+rateHHA(:,ib(j1))*hAnew(j1,idx);
                   else
                       tmp=rateHHA(:,ib(j))*hAnew(j,idx);
                   end
                   lambdaHH_new(:,idx)=lambdaHH(:,idx)+tmp;
                   lambdaHHA_new(:,idx)=lambdaHHA(:,idx)+tmp;
                   lambda_new(:,idx)=lambdaHH_new(ib,idx);
                   q=log(prob0(j,3)/sum(prob0(j,2:3))); %-log(probAIPnew(tRAnew(j)-tp));
               else % individual is already in asymptomatic infectiousness matrix as they recovered after t=0
                   actvAnew=actvA(actvA~=j); % remove from initially actively asymptomatic individuals
                   hAnew(j,idx1)=0; % remove previous infectiousness
%                    hAnew(j,idx2)=h4; % add new infectiousness
                   hAnew(j,idx2)=pold(6); % add new infectiousness
                   idx=union(idx1,idx2);
                   if ~mig_out || (mig_out && s<=rng(j,2)-1 && tRAnew(j)<=rng(j,2)-1)
                       tmp=rateHHA(:,ib(j))*(hAnew(j,idx)-hA(j,idx)); 
                   else % individual internally migrated and old or new recovery time is after migration time
                       if s<=rng(j,2)-1 && tRAnew(j)>rng(j,2)-1
                           RAobs2new=[RAobs2;j1];
                           idx4=rng(j,2)+1:min(min(tRAnew(j),rng(j1,2)),tmax);
                           hAnew(j1,idx4)=pold(6);
                           idx=[idx,idx4];
                           tmp=rateHHA(:,ib(j))*(hAnew(j,idx)-hA(j,idx))+rateHHA(:,ib(j1))*hAnew(j1,idx);
                       elseif s>rng(j,2)-1 && tRAnew(j)<=rng(j,2)-1
                           RAobs2actvAnew=RAobs2actvA(RAobs2actvA~=j1);
                           hAnew(j1,idx3)=0;
                           idx=[idx,idx3];
                           tmp=rateHHA(:,ib(j))*(hAnew(j,idx)-hA(j,idx))-rateHHA(:,ib(j1))*hA(j1,idx);
                       else %if s>rng(j,2)-1 && tRAnew(j)>rng(j,2)-1
                           RAobs2actvAnew=RAobs2actvA(RAobs2actvA~=j1);
                           RAobs2new=[RAobs2;j1];
                           hAnew(j1,idx3)=0; % remove infectiousness if currently still infected after moving
                           idx4=rng(j,2)+1:min(min(tRAnew(j),rng(j1,2)),tmax);
                           hAnew(j1,idx4)=pold(6);
                           idx=[idx,union(idx3,idx4)];
                           tmp=rateHHA(:,ib(j))*(hAnew(j,idx)-hA(j,idx))+rateHHA(:,ib(j1))*(hAnew(j1,idx)-hA(j1,idx));
                       end
                   end
                   lambdaHH_new(:,idx)=lambdaHH(:,idx)+tmp;
                   lambdaHHA_new(:,idx)=lambdaHHA(:,idx)+tmp;
                   lambda_new(:,idx)=lambdaHH_new(ib,idx);
                   q=log(prob0(j,2)/sum(prob0(j,2:3))); %+log(probAIP(s))-log(probAIPnew(tRAnew(j)-tp));
               end
               Susnew=[Sus;j]; %sort([Sus;j]); % add j to intially susceptible set
               S0new=[S0;j];
           end
%            q=q+log(numel(fwd))-log(numel(bck)); % calculate proposal ratio
%            q=q+log(probA(j,1)/sum(probA(j,bck+1)))-log(probA(j,tp+1)/sum(probA(j,fwd+1))); % calculate proposal ratio
%            q=q+log(probA(j,1))-log(probA(j,tp+1)); % calculate proposal ratio
           q=q+log(probA(j,1)/sum(probA(j,bck+1)))-log(probA(j,tp+1)/sum(probA(j,fwd+1))); % calculate proposal ratio
           
           Snew(j,1:min(tp-1,rng(j,2)-1))=1; % add susceptibility up to new infection time
           
           LL1new=L1(Snew,lambda_new);
%            LL2new=L2(lambda_new,pI,tEm);
%            LL4new=L4(lambda_new,pI,tAmnew);
           LL2new=L2(lambda_new,pold(7),tEm);
           LL4new=L4(lambda_new,pold(7),tAmnew);
%            LL5new=L5(age,S0new,actvAnew,prevAnew,lambda0,pI,p2);
           LL5new=L5(age,S0new,actvAnew,prevAnew,pold(5),pold(7),p2);
           LLnew=LL1new+LL2new+LL4new+LL5new;
           LLold=LL1old+LL2old+LL4old+LL5old;
           
           log_ap=LLnew-LLold+q; % calculate M-H acceptance probability
           
           if log_ap > log(rand)
               tA(j)=tp; tRA(j)=tRAnew(j); S(j,:)=Snew(j,:); tAm(j,:)=tAmnew(j,:);
               if mig_out && tp==tmax+1
                   S(j1,:)=Snew(j1,:);
                   tA(j1)=tmax+1; % change asymptomatic infection time for 2nd observation to tmax+1 for susceptibility
                   tRA(j1)=tmax+1; % change asymptomatic recovery time for 2nd observation to tmax+1 for susceptibility
               end
               hA=hAnew; lambda(:,idx)=lambda_new(:,idx); lambdaHH(:,idx)=lambdaHH_new(:,idx);
               lambdaHHA(:,idx)=lambdaHHA_new(:,idx);
               LL1old=LL1new; LL2old=LL2new;
               LL4old=LL4new; LL5old=LL5new;
               A1=A1new; nA1=nA1new; prevA=prevAnew; Sus=Susnew; actvA=actvAnew;
               S0=S0new;
               RAobs2actvA=RAobs2actvAnew; RAobs2=RAobs2new;
               acc_Arem=acc_Arem+1;
               acc_rem(j)=acc_rem(j)+1;
           else
               %tAnew(j)=tA(j); 
               tRAnew(j)=tRA(j); Snew(j,:)=S(j,:); tAmnew(j,:)=tAm(j,:);
               if mig_out && tp==tmax+1
                   Snew(j1,:)=S(j1,:);
               end
               hAnew=hA; lambda_new(:,idx)=lambda(:,idx); lambdaHH_new(:,idx)=lambdaHH(:,idx);
               lambdaHHA_new(:,idx)=lambdaHHA(:,idx);
               A1new=A1; nA1new=nA1; prevAnew=prevA; Susnew=Sus; actvAnew=actvA;
               S0new=S0;
               RAobs2actvAnew=RAobs2actvA; RAobs2new=RAobs2;
               rej_Arem=rej_Arem+1;
               rej_rem(j)=rej_rem(j)+1;
           end
       elseif t==tmax+1 % currently not asymptomatically infected before or during study (N.B. could be 1st or 2nd observation but only matters if it's 1st observation)
%           fwd=rng(j,1):rng(j,2)-1;
% %           tp=fwd(randi(numel(fwd)));
%           tp=fwd(randsample(numel(fwd),1,true,probA(j,fwd+1)));
%            fwd=0:tmax+1;
%            tp=fwd(sum(rand>=cprobA(j,:))+1);
           fwd=[rng(j,1):rng(j,2)-1,tmax+1];
           tp=fwd(sum(rand>=cumsum(probA(j,fwd+1))/sum(probA(j,fwd+1)))+1);
           if tp~=tmax+1 % calculate likelihood if proposed infection time is different
              if tp==0 % proposed asymptomatically infected before start of study (tp=0) (N.B. must be 1st observation)
%                   bck=[1:rng(j,2)-1,tmax+1]; % since rng(j,1) must be 0
%                   bck=0:tmax+1;
                  bck=[0:rng(j,2)-1,tmax+1];                  
                  if rand<prob0(j,3)/sum(prob0(j,2:3)) % recovered before t=0 
                      tRAnew(j)=0;
                      prevAnew=[prevA;j]; %sort([prevA;j]); % add j to previously asymptomatically infected set
                      idx=[];
                      q=-log(prob0(j,3)/sum(prob0(j,2:3)));
                  else % recovered after t=0
                      if mig_out
                          probAIPnew=[geopdf(0:rng(j1,2)-2,p2),1-geocdf(rng(j1,2)-2,p2)];
                          fwdRA=[1:rng(j1,2)-1,tmax+1];
                      else
                          probAIPnew=[geopdf(0:rng(j,2)-2,p2),1-geocdf(rng(j,2)-2,p2)];
                          fwdRA=[1:rng(j,2)-1,tmax+1];
                      end
%                       tRAnew(j)=fwdRA(randsample(numel(fwdRA),1,true,probAIPnew));
                      tRAnew(j)=fwdRA(sum(rand>=cumsum(probAIPnew))+1);
                      idx=1:min(min(tRAnew(j),rng(j,2)),tmax);
                      A1new=[A1;j];
                      nA1new=nA1+1;
                      actvAnew=[actvA;j];
%                       hAnew(j,idx)=h4;
                      hAnew(j,idx)=pold(6);
                      if mig_out && tRAnew(j)>rng(j,2)-1
                          RAobs2actvAnew=[RAobs2actvA;j1];
                          idx4=rng(j,2)+1:min(min(tRAnew(j),rng(j1,2)),tmax);
                          hAnew(j1,idx4)=pold(6);
                          idx=[idx,idx4];
                          tmp=rateHHA(:,ib(j))*hAnew(j,idx)+rateHHA(:,ib(j1))*hAnew(j1,idx);
                      else
                          tmp=rateHHA(:,ib(j))*hAnew(j,idx);
                      end                      
                      lambdaHH_new(:,idx)=lambdaHH(:,idx)+tmp;
                      lambdaHHA_new(:,idx)=lambdaHHA(:,idx)+tmp;
                      lambda_new(:,idx)=lambdaHH_new(ib,idx);
                      q=-log(prob0(j,2)/sum(prob0(j,2:3)));%+log(probAIPnew(tRAnew(j)));
                  end
                  Susnew=Sus(Sus~=j); % remove j from initially susceptible set
                  S0new=S0(S0~=j);
              else % proposed asymptomatically infected during study (tp in [1,tmax])
                  if mig_out
                      probAIPnew=[geopdf(0:rng(j1,2)-tp-2,p2),1-geocdf(rng(j1,2)-tp-2,p2)];
                      fwdRA=[tp+1:rng(j1,2)-1,tmax+1];
                  else
                      probAIPnew=[geopdf(0:rng(j,2)-tp-2,p2),1-geocdf(rng(j,2)-tp-2,p2)];
                      fwdRA=[tp+1:rng(j,2)-1,tmax+1];
                  end                  
%                   tRAnew(j)=fwdRA(randsample(numel(fwdRA),1,true,probAIPnew));
                  tRAnew(j)=fwdRA(sum(rand>=cumsum(probAIPnew))+1);
                  idx=tp+1:min(min(tRAnew(j),rng(j,2)),tmax);
    %               if rng(j,1)==0 % asymptomatic infection before study is possible as individual was alive/present
    %                   bck=[0,max(1,tp-M):tp-1,tp+1:min(rng(j,2),tp+M),tmax+1];
    %               else % individual was born or imigrated after start of study
    %                   bck=[max(rng(j,1),tp-M):tp-1,tp+1:min(rng(j,2),tp+M),tmax+1];
    %               end
%                   bck=[max(0,rng(j,1)):tp-1,tp+1:rng(j,2)-1,tmax+1];
                  if rng(j,1)==0 % asymptomatic infection before study is possible as individual was alive/present
                      bck=[0,max(1,tp-M):min(rng(j,2)-1,tp+M),tmax+1];
                  else % individual was born or imigrated after start of study
                      bck=[max(rng(j,1),tp-M):min(rng(j,2)-1,tp+M),tmax+1];
                  end
                  tAmnew(j,tp)=1; % add new infection time
                  A1new=[A1;j];
                  nA1new=nA1+1;
%                   hAnew(j,idx)=h4;
                  hAnew(j,idx)=pold(6);
                  if mig_out && tRAnew(j)>rng(j,2)-1
                      RAobs2new=[RAobs2;j1];
                      idx4=rng(j,2)+1:min(min(tRAnew(j),rng(j1,2)),tmax);
                      hAnew(j1,idx4)=pold(6);
                      idx=[idx,idx4];
                      tmp=rateHHA(:,ib(j))*hAnew(j,idx)+rateHHA(:,ib(j1))*hAnew(j1,idx);
                  else
                      tmp=rateHHA(:,ib(j))*hAnew(j,idx);
                  end
                  lambdaHH_new(:,idx)=lambdaHH(:,idx)+tmp;
                  lambdaHHA_new(:,idx)=lambdaHHA(:,idx)+tmp;
                  lambda_new(:,idx)=lambdaHH_new(ib,idx);
                  q=0; %-log(probAIPnew(tRAnew(j)-tp));
              end
    %           q=q+log(numel(fwd))-log(numel(bck)); % calculate proposal ratio
%               q=q+log(probA(j,tmax+2)/sum(probA(j,bck+1)))-log(probA(j,tp+1)/sum(probA(j,fwd+1))); % calculate proposal ratio
%               q=q+log(probA(j,tmax+2))-log(probA(j,tp+1)); % calculate proposal ratio
              q=q+log(probA(j,tmax+2)/sum(probA(j,bck+1)))-log(probA(j,tp+1)/sum(probA(j,fwd+1))); % calculate proposal ratio

              Snew(j,max(1,tp):rng(j,2)-1)=0; % remove susceptibility from new infection time onwards
              if mig_out
                  Snew(j1,rng(j,2):rng(j1,2)-1)=0; % remove susceptibility from 2nd observation now that individual is asymptomatically infected during 1st observation
              end

              LL1new=L1(Snew,lambda_new);
%               LL2new=L2(lambda_new,pI,tEm);
%               LL4new=L4(lambda_new,pI,tAmnew);
              LL2new=L2(lambda_new,pold(7),tEm);
              LL4new=L4(lambda_new,pold(7),tAmnew);
    %           LL5new=L5(age,S0new,actvAnew,prevAnew,lambda0,pI,p2);
              LL5new=L5(age,S0new,actvAnew,prevAnew,pold(5),pold(7),p2);
              LLnew=LL1new+LL2new+LL4new+LL5new;
              LLold=LL1old+LL2old+LL4old+LL5old;

              log_ap=LLnew-LLold+q; % calculate M-H acceptance probability

              if log_ap > log(rand)
                  tA(j)=tp; tRA(j)=tRAnew(j); S(j,:)=Snew(j,:); tAm(j,:)=tAmnew(j,:);
                  if mig_out
                      S(j1,:)=Snew(j1,:);
                      tA(j1)=tmax+2; % change asymptomatic infection time to dummy time for asymptomatic infection not possible now that individual is asymptomatically infected during 1st observation
                      tRA(j1)=tmax+2; % change asymptomatic recovery time to dummy time for asymptomatic infection not possible now that individual is asymptomatically infected during 1st observation
                  end
                  hA=hAnew; lambda(:,idx)=lambda_new(:,idx); lambdaHH(:,idx)=lambdaHH_new(:,idx);
                  lambdaHHA(:,idx)=lambdaHHA_new(:,idx);
                  LL1old=LL1new; LL2old=LL2new;
                  LL4old=LL4new; LL5old=LL5new;
                  A1=A1new; nA1=nA1new; prevA=prevAnew; Sus=Susnew; actvA=actvAnew;
                  S0=S0new;
                  RAobs2actvA=RAobs2actvAnew; RAobs2=RAobs2new;
                  acc_Aadd=acc_Aadd+1;
                  acc_add(j)=acc_add(j)+1;
              else
    %               tAnew(j)=tA(j);
                  tRAnew(j)=tRA(j); Snew(j,:)=S(j,:); tAmnew(j,:)=tAm(j,:);
                  if mig_out
                      Snew(j1,:)=S(j1,:);
                  end
                  hAnew=hA; lambda_new(:,idx)=lambda(:,idx); lambdaHH_new(:,idx)=lambdaHH(:,idx);
                  lambdaHHA_new(:,idx)=lambdaHHA(:,idx);
                  A1new=A1; nA1new=nA1; prevAnew=prevA; Susnew=Sus; actvAnew=actvA;
                  S0new=S0;
                  RAobs2actvAnew=RAobs2actvA; RAobs2new=RAobs2;
                  rej_Aadd=rej_Aadd+1;
                  rej_add(j)=rej_add(j)+1;
              end
           else % otherwise accept immediately as likelihood doesn't change (if tp=tmax+1)
               acc_Aadd=acc_Aadd+1;
           end
       else % currently asymptomatically infected between time 1 and tmax (N.B. could be 1st or 2nd observation and need to account for both)
%            if rng(j,1)==0
%                fwd=[0,max(1,t-M):t-1,t+1:min(rng(j,2),t+M),tmax+1];
%            else
%                fwd=[max(rng(j,1),t-M):t-1,t+1:min(rng(j,2),t+M),tmax+1];
%            end
%            fwd=[max(0,rng(j,1)):t-1,t+1:rng(j,2)-1,tmax+1];
% %            tp=fwd(randi(numel(fwd)));
%            tp=fwd(randsample(numel(fwd),1,true,probA(j,fwd+1)));
           if rng(j,1)==0
               fwd=[0,max(1,t-M):min(rng(j,2)-1,t+M),tmax+1];
           else
               fwd=[max(rng(j,1),t-M):min(rng(j,2)-1,t+M),tmax+1];
           end
           tp=fwd(sum(rand>=cumsum(probA(j,fwd+1))/sum(probA(j,fwd+1)))+1);
           tAmnew(j,t)=0; % remove old infection time
           hAnew(j,idx1)=0; % remove previous infectiousness
%            probAIP=[geopdf(0:tmax-t-1,p2),1-geocdf(tmax-t-1,p2)];
           if tp==0 % proposed asymptomatically infected before study (N.B. must be 1st observation)
%                bck=[1:rng(j,2)-1,tmax+1];
%                bck=0:tmax+1;
               bck=[0:rng(j,2)-1,tmax+1];
               Susnew=Sus(Sus~=j); % remove j from initially susceptible set
               S0new=S0(S0~=j);
               if rand<prob0(j,3)/sum(prob0(j,2:3)) % recovered before t=0
                   tRAnew(j)=0;
                   A1new=A1(A1~=j);
                   nA1new=nA1-1;
                   prevAnew=[prevA;j]; %sort([prevA;j]); % add j to previously asymptomatically infected set                   
                   idx=idx1; 
                   if mig_out && s>rng(j,2)-1
                       RAobs2new=RAobs2(RAobs2~=j1);
                       hAnew(j1,idx3)=0;
                       idx=[idx,idx3];
                       tmp=rateHHA(:,ib(j))*hA(j,idx)+rateHHA(:,ib(j1))*hA(j1,idx);
                   else
                       tmp=rateHHA(:,ib(j))*hA(j,idx);
                   end
                   lambdaHH_new(:,idx)=lambdaHH(:,idx)-tmp;
                   lambdaHHA_new(:,idx)=lambdaHHA(:,idx)-tmp;
                   lambda_new(:,idx)=lambdaHH_new(ib,idx);
                   q=-log(prob0(j,3)/sum(prob0(j,2:3)));%+log(probAIP(s-t))
               else % recovered after t=0
                   if mig_out
                       probAIPnew=[geopdf(0:rng(j1,2)-2,p2),1-geocdf(rng(j1,2)-2,p2)];
                       fwdRA=[1:rng(j1,2)-1,tmax+1];
                   else
                       probAIPnew=[geopdf(0:rng(j,2)-2,p2),1-geocdf(rng(j,2)-2,p2)];
                       fwdRA=[1:rng(j,2)-1,tmax+1];
                   end
%                    tRAnew(j)=fwdRA(randsample(numel(fwdRA),1,true,probAIPnew));
                   tRAnew(j)=fwdRA(sum(rand>=cumsum(probAIPnew))+1);
                   idx2=1:min(min(tRAnew(j),rng(j,2)),tmax);
                   actvAnew=[actvA;j];
%                    hAnew(j,idx2)=h4; % add new infectiousness
                   hAnew(j,idx2)=pold(6); % add new infectiousness
                   idx=union(idx1,idx2); % create union of indices for removed infectiousness and added infectiousness
                   if ~mig_out || (mig_out && s<=rng(j,2)-1 && tRAnew(j)<=rng(j,2)-1)
                       tmp=rateHHA(:,ib(j))*(hAnew(j,idx)-hA(j,idx));
                   else
                       if s<=rng(j,2)-1 && tRAnew(j)>rng(j,2)-1
                           RAobs2actvAnew=[RAobs2actvA;j1];
                           idx4=rng(j,2)+1:min(min(tRAnew(j),rng(j1,2)),tmax);
                           hAnew(j1,idx4)=pold(6);
                           idx=[idx,idx4];
                           tmp=rateHHA(:,ib(j))*(hAnew(j,idx)-hA(j,idx))+rateHHA(:,ib(j1))*hAnew(j1,idx);
                       elseif s>rng(j,2)-1 && tRAnew(j)<=rng(j,2)-1
                           RAobs2new=RAobs2(RAobs2~=j1);
                           hAnew(j1,idx3)=0;
                           idx=[idx,idx3];
                           tmp=rateHHA(:,ib(j))*(hAnew(j,idx)-hA(j,idx))-rateHHA(:,ib(j1))*hA(j1,idx);
                       else
                           RAobs2new=RAobs2(RAobs2~=j1);
                           RAobs2actvAnew=[RAobs2actvA;j1];
                           hAnew(j1,idx3)=0;
                           idx4=rng(j,2)+1:min(min(tRAnew(j),rng(j1,2)),tmax);
                           hAnew(j1,idx4)=pold(6);
                           idx=[idx,union(idx3,idx4)];
                           tmp=rateHHA(:,ib(j))*(hAnew(j,idx)-hA(j,idx))+rateHHA(:,ib(j1))*(hAnew(j1,idx)-hA(j1,idx));
                       end
                   end
                   lambdaHH_new(:,idx)=lambdaHH(:,idx)+tmp;
                   lambdaHHA_new(:,idx)=lambdaHHA(:,idx)+tmp;
                   lambda_new(:,idx)=lambdaHH_new(ib,idx);
                   q=-log(prob0(j,2)/sum(prob0(j,2:3)));%+log(probAIP(s-t))-log(probAIPnew(tRAnew(j)))
               end
           elseif tp==tmax+1 % proposed not asymptomatically infected before or during study
               tRAnew(j)=tmax+1;
%                bck=rng(j,1):rng(j,2)-1;
%                bck=0:tmax+1;
               bck=[rng(j,1):rng(j,2)-1,tmax+1];
               A1new=A1(A1~=j);
               nA1new=nA1-1;
               idx=idx1;
               if mig_out && s>rng(j,2)-1
                   RAobs2new=RAobs2(RAobs2~=j1);
                   hAnew(j1,idx3)=0;
                   idx=[idx,idx3];
                   tmp=rateHHA(:,ib(j))*hA(j,idx)+rateHHA(:,ib(j1))*hA(j1,idx);
               else
                   tmp=rateHHA(:,ib(j))*hA(j,idx);
               end
               lambdaHH_new(:,idx)=lambdaHH(:,idx)-tmp;
               lambdaHHA_new(:,idx)=lambdaHHA(:,idx)-tmp;
               lambda_new(:,idx)=lambdaHH_new(ib,idx);
               q=0; %log(probAIP(s-t));
               if mig_out
                   Snew(j1,rng(j,2):rng(j1,2)-1)=1; % add susceptibility to 2nd observation now that individual is not asymptomatically infected during 1st observation
               end
           else % proposed asymptomatically infected during study
               if mig_out
                   probAIPnew=[geopdf(0:rng(j1,2)-tp-2,p2),1-geocdf(rng(j1,2)-tp-2,p2)];
                   fwdRA=[tp+1:rng(j1,2)-1,tmax+1];
               else
                   probAIPnew=[geopdf(0:rng(j,2)-tp-2,p2),1-geocdf(rng(j,2)-tp-2,p2)];
                   fwdRA=[tp+1:rng(j,2)-1,tmax+1];
               end
%                tRAnew(j)=fwdRA(randsample(numel(fwdRA),1,true,probAIPnew));
               tRAnew(j)=fwdRA(sum(rand>=cumsum(probAIPnew))+1);
               idx2=tp+1:min(min(tRAnew(j),rng(j,2)),tmax);
%                if rng(j,1)==0 % asymptomatic infection before study is possible as individual was alive/present
%                    bck=[0,max(1,tp-M):tp-1,tp+1:min(rng(j,2),tp+M),tmax+1];
%                else % proposed asymptomatically infected during study and individual was born or imigrated after start of study
%                    bck=[max(rng(j,1),tp-M):tp-1,tp+1:min(rng(j,2),tp+M),tmax+1];
%                end
%                bck=[max(0,rng(j,1)):tp-1,tp+1:rng(j,2)-1,tmax+1];
               if rng(j,1)==0 % asymptomatic infection before study is possible as individual was alive/present
                   bck=[0,max(1,tp-M):min(rng(j,2)-1,tp+M),tmax+1];
               else % proposed asymptomatically infected during study and individual was born or imigrated after start of study
                   bck=[max(rng(j,1),tp-M):min(rng(j,2)-1,tp+M),tmax+1];
               end
               tAmnew(j,tp)=1; % add new infection time        
%                hAnew(j,idx2)=h4;
               hAnew(j,idx2)=pold(6);
               idx=union(idx1,idx2);
               if ~mig_out || (mig_out && s<=rng(j,2)-1 && tRAnew(j)<=rng(j,2)-1)
                   tmp=rateHHA(:,ib(j))*(hAnew(j,idx)-hA(j,idx));
               else
                   if s<=rng(j,2)-1 && tRAnew(j)>rng(j,2)-1
                       RAobs2new=[RAobs2;j1];
                       idx4=rng(j,2)+1:min(min(tRAnew(j),rng(j1,2)),tmax);
                       hAnew(j1,idx4)=pold(6);
                       idx=[idx,idx4];
                       tmp=rateHHA(:,ib(j))*(hAnew(j,idx)-hA(j,idx))+rateHHA(:,ib(j1))*hAnew(j1,idx);
                   elseif s>rng(j,2)-1 && tRAnew(j)<=rng(j,2)-1
                       RAobs2new=RAobs2(RAobs2~=j1);
                       hAnew(j1,idx3)=0;
                       idx=[idx,idx3];
                       tmp=rateHHA(:,ib(j))*(hAnew(j,idx)-hA(j,idx))-rateHHA(:,ib(j1))*hA(j1,idx);
                   else % if s>rng(j,2)-1 && tRAnew(j)>rng(j,2)-1
                       hAnew(j1,idx3)=0;
                       idx4=rng(j,2)+1:min(min(tRAnew(j),rng(j1,2)),tmax);
                       hAnew(j1,idx4)=pold(6);
                       idx=[idx,union(idx3,idx4)];
                       tmp=rateHHA(:,ib(j))*(hAnew(j,idx)-hA(j,idx))+rateHHA(:,ib(j1))*(hAnew(j1,idx)-hA(j1,idx));
                   end
               end
               lambdaHH_new(:,idx)=lambdaHH(:,idx)+tmp;
               lambdaHHA_new(:,idx)=lambdaHHA(:,idx)+tmp;
               lambda_new(:,idx)=lambdaHH_new(ib,idx);
               q=0; %log(probAIP(s-t))-log(probAIPnew(tRAnew(j)-tp));
           end
%            q=q+log(numel(fwd))-log(numel(bck)); % calculate proposal ratio
%            q=q+log(probA(j,t+1)/sum(probA(j,bck+1)))-log(probA(j,tp+1)/sum(probA(j,fwd+1))); % calculate proposal ratio
%            q=q+log(probA(j,t+1))-log(probA(j,tp+1)); % calculate proposal ratio
           q=q+log(probA(j,t+1)/sum(probA(j,bck+1)))-log(probA(j,tp+1)/sum(probA(j,fwd+1))); % calculate proposal ratio
           
           Snew(j,max(1,tp):t-1)=0; % remove susceptibility if new infection time is earlier
           Snew(j,t:min(tp-1,rng(j,2)-1))=1; % add susceptibility if new infection time is later
           
           LL1new=L1(Snew,lambda_new);
%            LL2new=L2(lambda_new,pI,tEm);
%            LL4new=L4(lambda_new,pI,tAmnew);
           LL2new=L2(lambda_new,pold(7),tEm);
           LL4new=L4(lambda_new,pold(7),tAmnew);
%            LL5new=L5(age,S0new,actvAnew,prevAnew,lambda0,pI,p2);
           LL5new=L5(age,S0new,actvAnew,prevAnew,pold(5),pold(7),p2);
           LLnew=LL1new+LL2new+LL4new+LL5new;
           LLold=LL1old+LL2old+LL4old+LL5old;
           
           log_ap=LLnew-LLold+q; % calculate M-H acceptance probability
           
           if log_ap > log(rand)
               tA(j)=tp; tRA(j)=tRAnew(j); S(j,:)=Snew(j,:); tAm(j,:)=tAmnew(j,:);
               if mig_out && tp==tmax+1
                   S(j1,:)=Snew(j1,:);
                   tA(j1)=tmax+1;
                   tRA(j1)=tmax+1;
               end
               hA=hAnew; lambda(:,idx)=lambda_new(:,idx); lambdaHH(:,idx)=lambdaHH_new(:,idx);
               lambdaHHA(:,idx)=lambdaHHA_new(:,idx);
               LL1old=LL1new; LL2old=LL2new;
               LL4old=LL4new; LL5old=LL5new;
               A1=A1new; nA1=nA1new; prevA=prevAnew; Sus=Susnew; actvA=actvAnew;
               S0=S0new;
               RAobs2actvA=RAobs2actvAnew; RAobs2=RAobs2new;
               acc_Amov=acc_Amov+1;
               acc_mov(j)=acc_mov(j)+1;
           else
%               tAnew(j)=tA(j);
               tRAnew(j)=tRA(j); Snew(j,:)=S(j,:); tAmnew(j,:)=tAm(j,:);
               if mig_out && tp==tmax+1
                   Snew(j1,:)=S(j1,:);
               end
               hAnew=hA; lambda_new(:,idx)=lambda(:,idx); lambdaHH_new(:,idx)=lambdaHH(:,idx);
               lambdaHHA_new(:,idx)=lambdaHHA(:,idx);
               A1new=A1; nA1new=nA1; prevAnew=prevA; Susnew=Sus; actvAnew=actvA;
               S0new=S0;
               RAobs2actvAnew=RAobs2actvA; RAobs2new=RAobs2;
               rej_Amov=rej_Amov+1;
               rej_mov(j)=rej_mov(j)+1;
           end
       end
       end
    end
%     % Update tAnew so code below works
%     tAnew=tA;
% 
%     %% UPDATE BETA AND ASYMPTOMATIC INFECTION TIMES JOINTLY
%     jmp=c1(k)*sqrt(ppvar1)*randn;
%     pnew(1)=pold(1)+jmp;
%     
%     if pnew(1)>plb(1) && pnew(1)<pub(1) % check if prior probability is non-zero before calculating log-likelihood
%         q=log(gampdf(pnew(1),prior_shape(1),prior_scale(1)))-log(gampdf(pold(1),prior_shape(1),prior_scale(1)));
%         
%         % NEEDS CORRECTING (constraining asx infctn times to be during
%         % study (>=1) is wrong)
%         A3=setdiff(A1,actvA);
%         tAnew(A3)=min(max(rng(A3,1)+1,round(tA(A3)+d(2)*jmp)),rng(A3,2));
%         tRAnew(A3)=min(max(rng(A3,1)+1,round(tRA(A3)+d(2)*jmp)),rng(A3,2));
%         j=find(tAnew==tmax+1 & tA<tmax+1);
%         A1new=A1(~ismember(A1,j));
%         nA1new=nA1-numel(j);
%         hAnew=false(n,tmax);
%         for i=1:nA1new
%             j=A1new(i);
%             hAnew(j,tAnew(j)+1:min(min(min(tRAnew(j),tEM(j)),tD(j)),tmax))=h4;
%         end
%         hAnew=sparse(hAnew);
%         tAmnew=false(n,tmax);
%         A2=find(tAnew>0 & tAnew<tmax+1);
%         tAmnew((tAnew(A2)-1)*n+A2)=1;
%         tAmnew=sparse(tAmnew);
%         
%         rateHHA_new=pnew(1)*KHH;
%         if ismember(4,u)
%             rateHHA_new=rateHHA_new+pold(4)*d0;
%         end
%         rateHH_new=rateHHA_new(:,IPNIA);
%         rateHHPA_new=rateHHA_new(:,PA);
%         lambdaHHA_new=rateHHA_new(:,A1new)*hAnew(A1new,:);
%         lambdaHHPA_new=rateHHPA_new*hPA;
%         lambdaHH_new=rateHH_new*h+lambdaHHPA_new+lambdaHHA_new+pold(3); % update whole infectious pressure
%         lambda_new=lambdaHH_new(ib,:);
%         
%         Snew(A3,:)=1-max(preB(A3,:),preIM(A3,:))-cumsum(tAm(A3,:),2);
%         
%         LL1new=L1(Snew,lambda_new);
%         LL2new=L2(lambda_new,pI,tEm);
%         LL4new=L4(lambda,pI,tAmnew);
% %         LL5new=L5(age,S0new,actvAnew,prevAnew,pold(5),pI,p2);
%         LLnew=LL1new+LL2new+LL4new;%+LL5new;
%         LLold=LL1old+LL2old+LL4old;%+LL5old;
%         
%         log_ap=LLnew-LLold+q; % calculate Metropolis-Hastings acceptance probability
%         
%         if log_ap > log(rand) % for acc_prob<=1, change if greater than rand
%             pold(1)=pnew(1); % keep updated parameter values
%             tA(A3)=tAnew(A3); tRA(A3)=tRAnew(A3); S(A3,:)=Snew(A3,:); tAm=tAmnew; hA=hAnew;
%             rateHHA=rateHHA_new; rateHH=rateHH_new;
%             lambda=lambda_new; lambdaHH=lambdaHH_new; % keep updated info
%             lambdaHHA=lambdaHHA_new; lambdaHHPA=lambdaHHPA_new;
%             A1=A1new; nA1=nA1new; 
%             LL1old=LL1new; LL2old=LL2new; LL4old=LL4new;
%             acc_beta=acc_beta+1;
%             c1(k+1)=(1+200/(1000+k))*c1(k);
%         else
%             pnew(1)=pold(1); % revert to old parameter value
%             tAnew(A3)=tA(A3); tRAnew(A3)=tRA(A3); Snew(A3,:)=S(A3,:); tAmnew=tAm; hAnew=hA;
%             lambda_new=lambda; lambdaHH_new=lambdaHH; % revert lambda_new to old value as only part of lambda_new is updated in infection time updates below
%             lambdaHHA_new=lambdaHHA; lambdaHHPA_new=lambdaHHPA;
%             A1new=A1; nA1new=nA1;
%             rej_beta=rej_beta+1;
%             c1(k+1)=(1+200/(1000+k))^(0.234/(0.234-1))*c1(k);
%         end
%     else % if prior probability is 0, reject immediately
%         pnew(1)=pold(1); % revert to old parameter value
%         rej_beta=rej_beta+1;
%         c1(k+1)=(1+200/(1000+k))^(0.234/(0.234-1))*c1(k);
%     end

    %% UPDATE ASYMPTOMATIC INFECTION AND RECOVERY TIMES OF PKDL CASES WITHOUT PRIOR KA
    for i=1:nPA
        j=PA(i);
        t=tA(j);
%         probDIP=[nbinpdf(0:tP(j)-1,r3,p3),1-nbincdf(tP(j)-1,r3,p3)];
%         tRAnew(j)=tP(j)-randsample(tP(j),1,true,probDIP);
        tRAnew(j)=tP(j)-nbinrnd(r3,p3);
        tp=tRAnew(j)-(geornd(p2)+1);
%         if tRAnew(j)<0
%             tRAnew(j)=0;
%         end
        
        if tp>=max(1,tB(j)+1) % asymptomatically infected after or during birth month
%             if tRA(j)>0 || tRAnew(j)>0
%                 if t>0 % currently asymptomatically infected during study
                    tAmnew(j,t)=0; % remove old infection time if it is after start of study
%                     if tp<=0 && tRAnew(j)>0 % proposed initially actively asymptomatically infected
%                         actvAPAnew=[actvAPA;j];
%                         S0PAnew=S0PA(S0PA~=j);
%                     elseif tRAnew(j)<=0 % proposed recovered before study
%                         prevAPAnew=[prevAPA;j];
%                         S0PAnew=S0PA(S0PA~=j);
%                     end
%                 else % currently asymptomatically infected before study (t<=0)
%                     if tRA(j)>0 % currently initially actively asymptomatically infected
%                         if tp>0 % proposed asymptomatically infected during study
%                             actvAPAnew=actvAPA(actvAPA~=j);
%                             S0PAnew=[S0PA;j];    
%                         elseif tRAnew(j)<=0 % proposed asymptomatically infected before study
%                             prevAPAnew=[prevAPA;j];
%                             actvAPAnew=actvAPA(actvAPA~=j);
%                         end
%                     else % tRA(j)<=0
%                         if tp<=0 && tRAnew(j)>0
%                             prevAPAnew=prevAPA(prevAPA~=j);
%                             actvAPAnew=[actvAPA;j];
%                         elseif tp>0
%                             prevAPAnew=prevAPA(prevAPA~=j);
%                             S0PAnew=[S0PA;j];
%                         end
%                     end
%                 end
%                 if tp>0
                    tAmnew(j,tp)=1; % add new infection time if it is after start of study
%                 end
                Snew(j,max(1,tp):t-1)=0; % remove old susceptible times if new infection time is earlier (from month indvdl is infctd up to but not incl. month before old infection time)
                Snew(j,max(1,t):tp-1)=1; % add new susceptible times if new infection time is later (up to but not incl. month indvdl is infctd)

                hPAnew(i,max(0,t)+1:tRA(j))=0; % remove old infectiousness
                hPAnew(i,max(0,tp)+1:tRAnew(j))=h40; % add new infectiousness
                erlrA=max(0,min(t,tp));
                ltrA=max(0,max(t,tp));
                erlrRA=max(0,min(tRA(j),tRAnew(j)));
                ltrRA=max(tRA(j),tRAnew(j));
                
                idx=[erlrA+1:ltrA,erlrRA+1:ltrRA];
                lambdaHHPA_new(:,idx)=lambdaHHPA(:,idx)-rateHHPA(:,i)*hPA(i,idx)+rateHHPA(:,i)*hPAnew(i,idx);
                lambdaHH_new(:,idx)=lambdaHH(:,idx)-rateHHPA(:,i)*hPA(i,idx)+rateHHPA(:,i)*hPAnew(i,idx);
                lambda_new(:,idx)=lambdaHH_new(ib,idx); % expand infectious pressure
                
                LL1new=L1(Snew,lambda_new);
%                 LL2new=L2(lambda_new,pI,tEm);
%                 LL4new=L4(lambda_new,pI,tAmnew);
                LL2new=L2(lambda_new,pold(7),tEm);
                LL4new=L4(lambda_new,pold(7),tAmnew);
%                 LL6new=L5(age,S0PAnew,actvAPAnew,prevAPAnew,lambda0,pI,p2);
                LL6new=L5(age,S0PAnew,actvAPAnew,prevAPAnew,pold(5),pold(7),p2);
                LLnew=LL1new+LL2new+LL4new+LL6new;
                LLold=LL1old+LL2old+LL4old+LL6old;
                
                log_ap=LLnew-LLold; % calculate Metropolis-Hastings acceptance probability
            
                if log_ap > log(rand)
                    tA(j)=tp; tRA(j)=tRAnew(j); S(j,:)=Snew(j,:); tAm(j,:)=tAmnew(j,:);
                    hPA=hPAnew; lambda(:,idx)=lambda_new(:,idx); lambdaHH(:,idx)=lambdaHH_new(:,idx);
                    lambdaHHPA(:,idx)=lambdaHHPA_new(:,idx);
                    LL1old=LL1new; LL2old=LL2new;
                    LL4old=LL4new; LL6old=LL6new;
                    prevAPA=prevAPAnew; S0PA=S0PAnew; actvAPA=actvAPAnew;
                    acc_PA=acc_PA+1;
                else
                    tRAnew(j)=tRA(j); Snew(j,:)=S(j,:); tAmnew(j,:)=tAm(j,:);
                    hPAnew=hPA; lambda_new(:,idx)=lambda(:,idx); lambdaHH_new(:,idx)=lambdaHH(:,idx);
                    lambdaHHPA_new(:,idx)=lambdaHHPA(:,idx);
                    prevAPAnew=prevAPA; S0PAnew=S0PA; actvAPAnew=actvAPA;
                    rej_PA=rej_PA+1;
                end
%             else % if both old and new recovery times are in/before start month then likelihood does not change so always accept
%                 tA(j)=tp; tRA(j)=tRAnew(j);
%                 acc_PA=acc_PA+1;
%             end
        else % reject immediately
            tRAnew(j)=tRA(j);
            rej_PA=rej_PA+1;
        end
    end

    %% NEED TO THINK ABOUT HOW TO CORRECT THE FOLLOWING TO ACCOUNT FOR IN-MIGRATION
    % E.G. SHOULD I RESTRICT MISSING ONSET TIMES TO BE AFTER IN-MIGRATION?
    %% UPDATE MISSING ONSET TIMES
    for i=1:nNO
        j=NO(i);
        tInew(j)=tR(j)-(nbinrnd(r0,p0)+1);
        
        if tInew(j)>tE(j) && tInew(j)>=tIlb(j) && tInew(j)<=tIub(j) % calculate log-likelihood if new onset time is after infection time, within infection time bounds, and before death
            tIj=tI(j);
            tIjnew=tInew(j);
            IPnew(j)=tIjnew-tE(j); % new incubation period
            m=(IPNIA==j);%ismember(IPNIA,j);
            hnew(m,tIj+1:tIjnew)=h0; % reduce infectiousness up to new onset time if it is later
            hnew(m,tIjnew+1:tIj)=1; % increase infectiousness from new onset time if it is earlier
            erlrI=min(tIj,tIjnew); % index of column for earlier onset time between old and new onset time
            ltrI=max(tIj,tIjnew); % index of column for later onset time between old and new onset time
            idx=erlrI+1:ltrI;
            lambdaHH_new(:,idx)=rateHH*hnew(:,idx)+lambdaHHPA(:,idx)+lambdaHHA(:,idx)+pold(3); % update infectious pressure
            lambda_new(:,idx)=lambdaHH_new(ib,idx); % expand infectious pressure

            LL1new=L1(S,lambda_new);
%             LL2new=L2(lambda_new,pI,tEm);
            LL2new=L2(lambda_new,pold(7),tEm);
            LL3new=L3(IPnew(I),r1,p1new);
%             LL4new=L4(lambda_new,pI,tAm);
            LL4new=L4(lambda_new,pold(7),tAm);
            LLnew=LL1new+LL2new+LL3new+LL4new;
            LLold=LL1old+LL2old+LL3old+LL4old;
            
            log_ap=LLnew-LLold; % calculate Metropolis-Hastings acceptance probability
            
            if log_ap > log(rand) % for acc_prob<=1, change if greater than rand
                IPold(j)=IPnew(j); tI(j)=tInew(j); h(m,:)=hnew(m,:); lambda(:,idx)=lambda_new(:,idx); lambdaHH(:,idx)=lambdaHH_new(:,idx);
                LL1old=LL1new; LL2old=LL2new; LL3old=LL3new; 
                LL4old=LL4new; % keep updated info
                acc_I=acc_I+1;
            else
                IPnew(j)=IPold(j); tInew(j)=tI(j); hnew(m,:)=h(m,:); lambda_new(:,idx)=lambda(:,idx); lambdaHH_new(:,idx)=lambdaHH(:,idx); % keep old values, don't change log-likelihood
                rej_I=rej_I+1;
            end
        else % otherwise reject immediately
            tInew(j)=tI(j);
            rej_I=rej_I+1;
        end
    end
 
    %% MOVE WHOLE INFECTION-TO-TREATMENT PERIOD OF CASES MISSING ONSET AND TREATMENT TIMES
    for i=1:nNONR
        j=NONR(i);
        tmp1=0;
        while tmp1==0
            tmp1=round(sqrt(ERvar)*randn);
        end
        tEnew(j)=tE(j)+tmp1;
        tInew(j)=tI(j)+tmp1;
        tRnew(j)=tR(j)+tmp1;
        
        % if KA onset and PKDL onset/death in same year and maximum symptom 
        % period currently chosen, or infection is at start of study/birth 
        % time and onset/recovery time are as late as they can be, reject 
        % immediately as infection-to-treatment period cannot be moved. N.B.
        % there is 1 KA case with simultaneous PKDL so PKDL onset can
        % happen before KA treatment, but we assume this is not the case in 
        % general, so enforce that PKDL onset occurs after KA onset
        if ~(tI(j)==tIlb(j) && tR(j)==min(min(tP(j)-1,tD(j)),tmax)) && ~(tE(j)==tB(j)+1 && (tI(j)==tIub(j) || tR(j)==min(min(tP(j)-1,tD(j)),tmax))) && tInew(j)>=tIlb(j) && tInew(j)<=tIub(j) && tEnew(j)>=tB(j)+1 && ~(tInew(j)>maxIP && tEnew(j)<1) && tRnew(j)<=min(min(tP(j)-1,tD(j)),tmax) % DO I NEED TO WORRY ABOUT EMIGRATION HERE?
            tEj=tE(j);
            tEjnew=tEnew(j);
            % calculate log-likelihood if an infection-treatment block move 
            % is possible; allow treatment month to be same as death month 
            % as cases were said to have died during treatment
            if tI(j)>maxIP && tInew(j)>maxIP % only update infection time and susceptibility matrices for cases with onset after initial window and after entry
                tEmnew(j,tEj)=0; % remove old infection time
                tEmnew(j,tEjnew)=1; % add new infection time
                Snew(j,tEjnew:tEj-1)=0; % remove old susceptible times if new infection time is earlier (from month indvdl is infctd up to but not incl. month before old infection time)
%                 Snew(j,tEj:tEjnew-1)=sigma(j,tEj:tEjnew-1); % add new susceptible times if new infection time is later (up to but not incl. month indvdl is infctd)
                Snew(j,tEj:tEjnew-1)=1; % add new susceptible times if new infection time is later (up to but not incl. month indvdl is infctd)
            elseif tI(j)>maxIP && tInew(j)<=maxIP
                tEmnew(j,tEj)=0;
                Snew(j,:)=0;
            elseif tI(j)<=maxIP && tInew(j)>maxIP
                tEmnew(j,tEjnew)=1;
                Snew(j,1:tEjnew)=1;
            end
            
            tIj=tI(j);
            tIjnew=tInew(j);
            tRj=tR(j);
            tRjnew=tRnew(j);            
            m=(IPNIA==j);%ismember(IPNIA,j);
            hnew(m,max(0,tEj)+1:tIj)=0; % remove old presymptomatic infectiousness
            hnew(m,tIj+1:tRj)=0; % remove old symptomatic infectiousness
            hnew(m,max(0,tEjnew)+1:tIjnew)=h0; % add new infectiousness
            hnew(m,tIjnew+1:tRjnew)=1; % increase infectiousness from new onset time if it is earlier
  
            erlrE=max(0,min(tEj,tEjnew)); % index of column for earlier infctn time between old and new infection time
            ltrE=max(0,max(tEj,tEjnew)); % index of column for later infctn time between old and new infection time
            erlrI=min(tIj,tIjnew); % index of column for earlier onset time between old and new onset time
            ltrI=max(tIj,tIjnew); % index of column for later onset time between old and new onset time
            erlrR=min(tRj,tRjnew); % index of column for earlier treatment time between old and new treatment time
            ltrR=max(tRj,tRjnew); % index of column for later treatment time between old and new treatment time
            
            idx=[erlrE+1:ltrE,erlrI+1:ltrI,erlrR+1:ltrR];
            lambdaHH_new(:,idx)=rateHH*hnew(:,idx)+lambdaHHPA(:,idx)+lambdaHHA(:,idx)+pold(3); % update infectious pressure
            lambda_new(:,idx)=lambdaHH_new(ib,idx); % expand infectious pressure
            
            LL1new=L1(Snew,lambda_new);
%             LL2new=L2(lambda_new,pI,tEmnew);
%             LL4new=L4(lambda_new,pI,tAm);
            LL2new=L2(lambda_new,pold(7),tEmnew);
            LL4new=L4(lambda_new,pold(7),tAm);
            LLnew=LL1new+LL2new+LL4new;
            LLold=LL1old+LL2old+LL4old;
            
            log_ap=LLnew-LLold; % calculate Metropolis-Hastings acceptance probability
            
            if log_ap > log(rand) % for acc_prob<=1, change if greater than rand
                tE(j)=tEnew(j); S(j,:)=Snew(j,:); tEm(j,:)=tEmnew(j,:); tI(j)=tInew(j); tR(j)=tRnew(j);
                h(m,:)=hnew(m,:); lambda(:,idx)=lambda_new(:,idx); lambdaHH(:,idx)=lambdaHH_new(:,idx);
                LL1old=LL1new; LL2old=LL2new; 
                LL4old=LL4new; % keep updated info
                acc_ERmove=acc_ERmove+1;
            else
                tEnew(j)=tE(j); Snew(j,:)=S(j,:); tEmnew(j,:)=tEm(j,:); tInew(j)=tI(j); tRnew(j)=tR(j);
                hnew(m,:)=h(m,:); lambda_new(:,idx)=lambda(:,idx); lambdaHH_new(:,idx)=lambdaHH(:,idx); % keep old values, don't change log-likelihood
                rej_ERmove=rej_ERmove+1;
            end
        else % otherwise reject immediately
            tEnew(j)=tE(j); tInew(j)=tI(j); tRnew(j)=tR(j);
            rej_ERmove=rej_ERmove+1;
        end
    end
    
    %% UPDATE MISSING TREATMENT TIMES
    for i=1:nONR
        j=ONR(i);
        tRnew(j)=tI(j)+nbinrnd(r0,p0)+1;
              
        if tRnew(j)<=min(min(tP(j)-1,tD(j)),tmax) % calculate log-likelihood if new recovery time is before PKDL onset/death/end of study
            tRj=tR(j);
            tRjnew=tRnew(j);
            m=(IPNIA==j);%ismember(IPNIA,j);
            hnew(m,tRj+1:tRjnew)=1; % increase infectiousness up to new treatment time if it is later
            hnew(m,tRjnew+1:tRj)=0; % reduce infectiousness from new treatment time if it is earlier
            erlrR=min(tRj,tRjnew); % index of column for earlier treatment time between old and new treatment time
            ltrR=max(tRj,tRjnew); % index of column for later treatment time between old and new treatment time
            idx=erlrR+1:ltrR;
            lambdaHH_new(:,idx)=rateHH*hnew(:,idx)+lambdaHHPA(:,idx)+lambdaHHA(:,idx)+pold(3); % update infectious pressure
            lambda_new(:,idx)=lambdaHH_new(ib,idx); % expand infectious pressure
            
            LL1new=L1(S,lambda_new);
%             LL2new=L2(lambda_new,pI,tEm);
%             LL4new=L4(lambda_new,pI,tAm);
            LL2new=L2(lambda_new,pold(7),tEm);
            LL4new=L4(lambda_new,pold(7),tAm);
            LLnew=LL1new+LL2new+LL4new;
            LLold=LL1old+LL2old+LL4old;
            
            log_ap=LLnew-LLold; % calculate Metropolis-Hastings acceptance value
            
            if log_ap > log(rand) % for acc_prob<=1, change if greater than rand
                tR(j)=tRnew(j); h(m,:)=hnew(m,:); lambda(:,idx)=lambda_new(:,idx); lambdaHH(:,idx)=lambdaHH_new(:,idx);
                LL1old=LL1new; LL2old=LL2new; 
                LL4old=LL4new; % keep updated info
                acc_R=acc_R+1;
            else
                tRnew(j)=tR(j); hnew(m,:)=h(m,:); lambda_new(:,idx)=lambda(:,idx); lambdaHH(:,idx)=lambdaHH_new(:,idx); % keep old values, don't change log-likelihood
                rej_R=rej_R+1;
            end
        else % otherwise reject immediately
            tRnew(j)=tR(j);
            rej_R=rej_R+1;
        end
    end
    
    %% MOVE WHOLE ONSET-TO-TREATMENT PERIOD OF CASES WITH ACTIVE KA AT START OF STUDY MISSING ONSET AND TREATMENT TIMES
    for i=1:nANONR
        j=ANONR(i);
        tmp2=0;
        while tmp2==0
            tmp2=round(sqrt(ERvar)*randn);
        end
        tInew(j)=tI(j)+tmp2;
        tRnew(j)=tInew(j)+nbinrnd(r0,p0)+1; % draw new onset-to-treatment period from onset-to-treatment distribution
            
        % if KA onset and PKDL onset/death in same year and maximum symptom 
        % period currently chosen, reject immediately as 
        % infection-to-treatment period cannot be moved
        if ~(tI(j)==tIlbA(j) && tR(j)==min(min(tP(j)-1,tD(j)),tmax)) && tInew(j)>=tIlbA(j) && tInew(j)<=tIubA(j) && tRnew(j)<=min(min(tP(j)-1,tD(j)),tmax)
            % calculate log-likelihood if an infection-treatment block move 
            % is possible; allow treatment month to be same as death month 
            % as cases were said to have died during treatment                       
            tRj=tR(j);
            tRjnew=tRnew(j);
            
            if tRj>0 || tRjnew>0 % update likelihood if old or new recovery time is after start month
                m=(IPNIA==j);%ismember(IPNIA,j);
                hnew(m,max(0,tRj)+1:tRjnew)=1; % increase infectiousness up to new treatment time if it is later
                hnew(m,max(0,tRjnew)+1:tRj)=0; % reduce infectiousness from new treatment time if it is later
                erlrR=max(0,min(tRj,tRjnew)); % index of column for earlier treatment time between old and new treatment time
                ltrR=max(tRj,tRjnew); % index of column for later treatment time between old and new treatment time
                
                idx=erlrR+1:ltrR;
                lambdaHH_new(:,idx)=rateHH*hnew(:,idx)+lambdaHHPA(:,idx)+lambdaHHA(:,idx)+pold(3); % update infectious pressure
                lambda_new(:,idx)=lambdaHH_new(ib,idx); % expand infectious pressure
                
                LL1new=L1(S,lambda_new);
%                 LL2new=L2(lambda_new,pI,tEm);
%                 LL4new=L4(lambda_new,pI,tAm);
                LL2new=L2(lambda_new,pold(7),tEm);
                LL4new=L4(lambda_new,pold(7),tAm);
                LLnew=LL1new+LL2new+LL4new;
                LLold=LL1old+LL2old+LL4old;
                
                log_ap=LLnew-LLold; % calculate Metropolis-Hastings acceptance probability
                
                if log_ap > log(rand) % for acc_prob<=1, change if greater than rand
                    tI(j)=tInew(j); tR(j)=tRnew(j);
                    h(m,:)=hnew(m,:); lambda(:,idx)=lambda_new(:,idx); lambdaHH(:,idx)=lambdaHH_new(:,idx); 
                    LL1old=LL1new; LL2old=LL2new; 
                    LL4old=LL4new; % keep updated info
                    acc_AIRmove=acc_AIRmove+1;
                else
                    tInew(j)=tI(j); tRnew(j)=tR(j);
                    hnew(m,:)=h(m,:); lambda_new(:,idx)=lambda(:,idx); lambdaHH_new(:,idx)=lambdaHH(:,idx); % keep old values, don't change log-likelihood
                    rej_AIRmove=rej_AIRmove+1;
                end
            else % if both old and new recovery times are in/before start month then likelihood does not change so always accept
                tI(j)=tInew(j); tR(j)=tRnew(j);
                acc_AIRmove=acc_AIRmove+1;
            end
        else
            tInew(j)=tI(j); tRnew(j)=tR(j);
            rej_AIRmove=rej_AIRmove+1;
        end
    end

    %% UPDATE MISSING TREATMENT TIMES OF CASES WITH ACTIVE KA AT START OF STUDY
    for i=1:nAONR
        j=AONR(i);
        tRnew(j)=tI(j)+nbinrnd(r0,p0)+1;
        
        if tRnew(j)<=min(min(tP(j)-1,tD(j)),tmax) % calculate log-likelihood if new recovery time is before PKDL onset/death/end of study
            tRj=tR(j);
            tRjnew=tRnew(j);           
            if tRj>0 || tRjnew>0 % calculate log-likelihood if old or new recovery time is after start month
                m=(IPNIA==j);%ismember(IPNIA,j);
                hnew(m,max(0,tRj)+1:tRjnew)=1; % add infectiousness up to new treatment time if it is later
                hnew(m,max(0,tRjnew)+1:tRj)=0; % remove infectiousness from new treatment time if it is earlier
                erlrR=max(0,min(tRj,tRjnew)); % index of column for earlier treatment time between old and new treatment time
                ltrR=max(tRj,tRjnew); % index of column for later treatment time between old and new treatment time
                idx=erlrR+1:ltrR;
                lambdaHH_new(:,idx)=rateHH*hnew(:,idx)+lambdaHHPA(:,idx)+lambdaHHA(:,idx)+pold(3); % update infectious pressure
                lambda_new(:,idx)=lambdaHH_new(ib,idx); % expand infectious pressure                    
                
                LL1new=L1(S,lambda_new);
%                 LL2new=L2(lambda_new,pI,tEm);
%                 LL4new=L4(lambda_new,pI,tAm);
                LL2new=L2(lambda_new,pold(7),tEm);
                LL4new=L4(lambda_new,pold(7),tAm);
                LLnew=LL1new+LL2new+LL4new;
                LLold=LL1old+LL2old+LL4old;                
                
                log_ap=LLnew-LLold; % calculate Metropolis-Hastings acceptance value
                
                if log_ap > log(rand) % for acc_prob<=1, change if greater than rand
                    tR(j)=tRnew(j); h(m,:)=hnew(m,:); lambda(:,idx)=lambda_new(:,idx); lambdaHH(:,idx)=lambdaHH_new(:,idx); 
                    LL1old=LL1new; LL2old=LL2new;
                    LL4old=LL4new; % keep updated info
                    acc_AR=acc_AR+1;
                else
                    tRnew(j)=tR(j); hnew(m,:)=h(m,:); lambda_new(:,idx)=lambda(:,idx); lambdaHH_new(:,idx)=lambdaHH(:,idx); % keep old values, don't change log-likelihood
                    rej_AR=rej_AR+1;
                end
            else % if both old and new recovery times are in/before start month then likelihood does not change so always accept
                tR(j)=tRnew(j);
                acc_AR=acc_AR+1;
            end
        else % otherwise reject immediately
            tRnew(j)=tR(j);
            rej_AR=rej_AR+1;
        end
    end
    
    %% MOVE WHOLE RELAPSE PERIOD OF RELAPSE CASES MISSING RELAPSE AND RELAPSE TREATMENT TIMES
    for i=1:nRLNO
       j=RLNO(i);
       tRLnew(j)=tR(j)+geornd(p4)+1;
       tRLRnew(j)=tRLnew(j)+nbinrnd(r0,p0)+1;
       
       if tRLnew(j)<=min(min(min(tEM(j),tP(j)),tD(j))-2,tmax) && tRLRnew(j)<=min(min(min(tEM(j),tP(j)),tD(j))-1,tmax)
           tRLj=tRL(j);
           tRLjnew=tRLnew(j);
           tRLRj=tRLR(j);
           tRLRjnew=tRLRnew(j);
           m=(IPNIA==j);
           hnew(m,tRLj+1:tRLRj)=0; % remove current infectiousness
           hnew(m,tRLjnew+1:tRLRjnew)=1; % add new infectiousness
           erlrRL=min(tRLj,tRLjnew); % index of column for earlier onset time between old and new onset time
           ltrRL=max(tRLj,tRLjnew); % index of column for later onset time between old and new onset time
           erlrRLR=min(tRLRj,tRLRjnew);
           ltrRLR=max(tRLRj,tRLRjnew);
           
           idx=[erlrRL+1:ltrRL,erlrRLR+1:ltrRLR];
           lambdaHH_new(:,idx)=rateHH*hnew(:,idx)+lambdaHHPA(:,idx)+lambdaHHA(:,idx)+pold(3); % update infectious pressure
           lambda_new(:,idx)=lambdaHH_new(ib,idx); % expand infectious pressure
           
           LL1new=L1(S,lambda_new);
%            LL2new=L2(lambda_new,pI,tEm);
           LL2new=L2(lambda_new,pold(7),tEm);
%            LL4new=L4(lambda_new,pI,tAm);
           LL4new=L4(lambda_new,pold(7),tAm);
           LLnew=LL1new+LL2new+LL4new;
           LLold=LL1old+LL2old+LL4old;
           
           log_ap=LLnew-LLold; % calculate Metropolis-Hastings acceptance probability
           
           if log_ap > log(rand) % for acc_prob<=1, change if greater than rand
               tRL(j)=tRLnew(j); tRLR(j)=tRLRnew(j); h(m,:)=hnew(m,:); lambda(:,idx)=lambda_new(:,idx); lambdaHH(:,idx)=lambdaHH_new(:,idx);
               LL1old=LL1new; LL2old=LL2new;
               LL4old=LL4new; % keep updated info
               acc_RLNO=acc_RLNO+1;
           else
               tRLnew(j)=tRL(j); tRLRnew(j)=tRLR(j); hnew(m,:)=h(m,:); lambda_new(:,idx)=lambda(:,idx); lambdaHH_new(:,idx)=lambdaHH(:,idx); % keep old values, don't change log-likelihood
               rej_RLNO=rej_RLNO+1;
           end
       else % otherwise reject immediately
           tRLnew(j)=tRL(j); tRLRnew(j)=tRLR(j);
           rej_RLNO=rej_RLNO+1;
       end
    end
    %% UPDATE MISSING RELAPSE TREATMENT TIMES
    for i=1:nRLO
       j=RLO(i);
       tRLRnew(j)=tRL(j)+nbinrnd(r0,p0)+1;
       
       if tRLRnew(j)<=min(min(min(tEM(j),tP(j)),tD(j))-1,tmax)
           tRLRj=tRLR(j);
           tRLRjnew=tRLRnew(j);
           m=(IPNIA==j);
           hnew(m,tRLRj+1:tRLRjnew)=1; % add infectiousness if new relapse treatment time is later
           hnew(m,tRLRjnew+1:tRLRj)=0; % remove infectiousness if new relapse treatment time is earlier
           erlrRLR=min(tRLRj,tRLRjnew); % index of column for earlier onset time between old and new onset time
           ltrRLR=max(tRLRj,tRLRjnew); % index of column for later onset time between old and new onset time
           idx=erlrRLR+1:ltrRLR;
           lambdaHH_new(:,idx)=rateHH*hnew(:,idx)+lambdaHHPA(:,idx)+lambdaHHA(:,idx)+pold(3); % update infectious pressure
           lambda_new(:,idx)=lambdaHH_new(ib,idx); % expand infectious pressure
           
           LL1new=L1(S,lambda_new);
%            LL2new=L2(lambda_new,pI,tEm);
           LL2new=L2(lambda_new,pold(7),tEm);
%            LL4new=L4(lambda_new,pI,tAm);
           LL4new=L4(lambda_new,pold(7),tAm);
           LLnew=LL1new+LL2new+LL4new;
           LLold=LL1old+LL2old+LL4old;
           
           log_ap=LLnew-LLold; % calculate Metropolis-Hastings acceptance probability
           
           if log_ap > log(rand) % for acc_prob<=1, change if greater than rand
               tRLR(j)=tRLRnew(j); h(m,:)=hnew(m,:); lambda(:,idx)=lambda_new(:,idx); lambdaHH(:,idx)=lambdaHH_new(:,idx);
               LL1old=LL1new; LL2old=LL2new;
               LL4old=LL4new; % keep updated info
               acc_RLO=acc_RLO+1;
           else
               tRLRnew(j)=tRLR(j); hnew(m,:)=h(m,:); lambda_new(:,idx)=lambda(:,idx); lambdaHH_new(:,idx)=lambdaHH(:,idx); % keep old values, don't change log-likelihood
               rej_RLO=rej_RLO+1;
           end
       else % otherwise reject immediately
           tRLRnew(j)=tRLR(j);
           rej_RLO=rej_RLO+1;
       end
    end
    %% SAVE PARAMETER VALUES AND PLOT PROGRESS
    % Save parameters, log-likelihood, and asymptomatic infection periods
%     p(k,:)=pold;
    p(k+1,:)=pold;
    p1(k)=p1new;
    K0(k)=K0old;
%     LL(k)=LLold;
    terms(k,:)=[LL1old,LL2old,LL3old,LL4old,LL5old,LL6old];
%     terms(k,:)=[LL1old,LL2old,LL3old,LL4old,LL5old];
    LL(k)=sum(terms(k,:));
    IPs(:,k)=IPold(I);
    tEs(:,k)=tE(I);
%     tAs(:,k)=tA;
%     tRAs(:,k)=tRA;
    tAs([SusA;PA],k)=tA([SusA;PA]);
    tRAs([SusA;PA],k)=tRA([SusA;PA]);
    tIsNONR(:,k)=tI(NONR);
    tRsNONR(:,k)=tR(NONR);
    tIsRNO(:,k)=tI(RNO);
    tRsONR(:,k)=tR(ONR);
    tIsANONR(:,k)=tI(ANONR);
    tRsANONR(:,k)=tR(ANONR);
    tRsAONR(:,k)=tR(AONR);
    tRLsRLO(:,k)=tRL(RLO);
    tRLRsRLO(:,k)=tRLR(RLO);
    tRLsRLNO(:,k)=tRL(RLNO);
    tRLRsRLNO(:,k)=tRLR(RLNO);
    
    % Update empirical mean and covariance for proposal distribution
%     [ppmean,ppvar]=updateMeanAndCov(ppmean,ppvar,pold,k);
%     [ppmean,ppvar]=updateMeanAndCovSpencer(ppmean,ppvar,p,k,burnin);
%     [ppmean,ppvar]=updateMeanAndCovSpencer(ppmean,ppvar,p,k,100);
    [ppmean,ppvar]=updateMeanAndCovSpencer(ppmean,ppvar,p,k,1000);
    lambda_mean=updateMeanInfctsPressure(lambda_mean,lambda,k);
    
    % Calculate acceptance rates for parameters and infection time moves
    acc_rate_p=acc_p/(acc_p+rej_p);
    acc_rate_E=acc_E/(acc_E+rej_E);
    acc_rate_I=acc_I/(acc_I+rej_I);
    acc_rate_ERmove=acc_ERmove/(acc_ERmove+rej_ERmove);
    acc_rate_R=acc_R/(acc_R+rej_R);
    acc_rate_AIRmove=acc_AIRmove/(acc_AIRmove+rej_AIRmove);
    acc_rate_AR=acc_AR/(acc_AR+rej_AR);
    acc_rate_RLNO=acc_RLNO/(acc_RLNO+rej_RLNO);
    acc_rate_RLO=acc_RLO/(acc_RLO+rej_RLO);
    acc_rate_Aadd=acc_Aadd/(acc_Aadd+rej_Aadd);
    acc_rate_Arem=acc_Arem/(acc_Arem+rej_Arem);
    acc_rate_Amov=acc_Amov/(acc_Amov+rej_Amov);
    acc_rate_PA=acc_PA/(acc_PA+rej_PA);
%     acc_rate_beta=acc_beta/(acc_beta+rej_beta);

%     lambdaHH1=rateHH*h+rateHHA(:,A1)*hA(A1,:)+rateHHA(:,PA)*hPA+pold(3);
%     error_lambda(k)=sum(sum(abs(lambdaHH1-lambdaHH)));
%     error_lambda(k)

    % Print output
    if mod(k,1e3)==0 || k==niters % every 1000 iterations
        fprintf('Iteration %d done.\n', k); % display iteration number
        fprintf('Current likelihood=%6.4g\n', LL(k));
        fprintf('LL1=%6.4g\n', terms(k,1));
        fprintf('LL2=%6.4g\n', terms(k,2));
        fprintf('LL3=%6.4g\n', terms(k,3));
        fprintf('LL4=%6.4g\n', terms(k,4));
        fprintf('LL5=%6.4g\n', terms(k,5));
%         fprintf('LL6=%6.4g\n', terms(k,6));
        fprintf('Parameter block update acceptance rate=%5.3f%%\n',100*acc_rate_p);
        
%         % Plot output
%         if plotOutpt && k>burnin % ignore burn-in period
%             z=burnin+1:k; % iterations to plot
%             figure(4);
%             [mode_p,HPDI,mode_p1,HPDI1]=PlotOutput(z,LL,p,np,pname,prior_mean,p1,a,b,n,tmax,I,RpreD,DpreR,tI,tR,tD,tRLm,tRLRm,nbins,scrnsz);            
%             figure(5);
%             PlotTrace(z,p,np,pname,p1,mode_p,HPDI,mode_p1,HPDI1,scrnsz)
%             drawnow
%             for j=1:np
%                 fprintf(['mode ' pname{j} '=%6.4g\n'], mode_p(j));
%             end
%             fprintf('mode p1=%6.4g\n', mode_p1);
%         end
    end
end
 
%% Display acceptance rates
fprintf('Parameter block update acceptance rate is %5.3f%%.\n',100*acc_rate_p);
fprintf('Infection time move acceptance rate is %5.3f%%.\n',100*acc_rate_E);
fprintf('Onset time move acceptance rate is %5.3f%%.\n',100*acc_rate_I);
fprintf('Infection period move acceptance rate is %5.3f%%.\n',100*acc_rate_ERmove);
fprintf('Treatment time move acceptance rate is %5.3f%%.\n',100*acc_rate_R);
fprintf('Active KA case onset-treatment period move acceptance rate is %5.3f%%.\n',100*acc_rate_AIRmove);
fprintf('Active KA case treatment time move acceptance rate is %5.3f%%.\n',100*acc_rate_AR);
fprintf('Relapse KA case relapse period move acceptance rate is %5.3f%%.\n',100*acc_rate_RLNO);
fprintf('Relapse KA case treatment time move acceptance rate is %5.3f%%.\n',100*acc_rate_RLO);
fprintf('Asymptomatic infection removal acceptance rate is %5.3f%%.\n',100*acc_rate_Arem);
fprintf('Asymptomatic infection time move acceptance rate is %5.3f%%.\n',100*acc_rate_Amov);
fprintf('Asymptomatic infection addition acceptance rate is %5.3f%%.\n',100*acc_rate_Aadd);
fprintf('Asymptomatic infection period move for PKDL cases w/o VL acceptance rate is %5.3f%%.\n',100*acc_rate_PA);
% fprintf('Joint beta asymptomatic infection time update acceptance rate is %5.3f%%.\n',100*acc_rate_beta);

% Remove first row (initial values) of p
p=p(2:end,:);

clear data dHH dHHsqrd KHH KHH_new rateHHA rateHHA_new rateHH rateHH_new
% save(rslts)
save(rslts,'-v7.3')
