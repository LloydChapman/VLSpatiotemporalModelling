function VLStmMCMC2010FastLklhdMgrtn(data,r1,p10,a,b,p2,u,beta0,alpha0,epsilon0,delta0,lambda0,h0,h1,h2,h3,hmssng,h40,pI0,typ,niters,plotOutpt,rslts,para)

%% Set default input parameters if no inputs are supplied
if nargin==0
    load('C:\Users\timpo\OneDrive - University of Warwick\taubern_baybern\raw_data_plus_cleaning\matlab_bayesianmodel\data_final2.mat') % database
    % Parameters for negative binomial incubation period distribution
    r1=3; % shape parameter
    mu=5; % mean incubation period in months
    p10=r1/(mu-1+r1); % initial value for 'success' parameter
    % Shape parameters for beta prior for p1
    a=22; 
    b=28;
    p2=1/5; % geometric Geom(p2) asymptomatic infection period distribution parameter
    u=1:4; % transmission parameters to estimate
    beta0=3; % spatial transmission rate constant
    alpha0=100; % distance scale factor for spatial kernel;
    epsilon0=1e-3; % background transmission rate
    delta0=1e-2; % additional within-HH transmission rate
    % Estimate historical asymptomatic infection rate by fitting initial status model to LST data from 2002
    s1=load('C:\Users\timpo\OneDrive - University of Warwick\taubern_baybern\1999_2004dataset\data_final.mat');
    [pars,~]=FitCatModLST3(s1.data,p2);
    lambda0=pars/12; % historical asymptomatic infection rate
    h0=0.02; % relative infectiousness of pre-symptomatics
    % Relative infectiousness of different forms of PKDL
    h1=9/26/(10/15); % macular/papular
    h3=18/21/(10/15); % nodular
    h2=(h1+h3)/2; % plaque (assumed halfway between macular/papular and nodular)
    h40=0.02; % relative infectiousness of asymptomatic individuals
    pI0=0.15; % proportion of infections that lead to VL
    hmssng=(101/138*h1+31/138*h2+6/138*h3); % unexamined - use average relative infectiousness of examined PKDL cases
    typ='Exp'; % type of transmission kernel
    niters=1e5; % no. of MCMC iterations
    plotOutpt=true; % plot output (true) or not (false)
    rslts='MCMC_NBIP_PKDL_ASX.mat'; % name of results file for output
    para=1; % paras to include in the data
end
 
%% LOAD DATA 
% Select data for para
data=data(ismember(data.PARA,para),:);

% Rename longitude and latitude variables
data.Properties.VariableNames{'HHNEWLNG'}='longitude';
data.Properties.VariableNames{'HHNEWLAT'}='latitude';
 
% Number of individuals
n=size(data,1);

% Assume cases stop being infectious shortly after commencing treatment
durRX=0; % months
 
% Set start year and start month
startyr=2002;
startmo=1;
 
% Set end year and month
endyr=2010;
endmo=12;
 
% Set maximum incubation period for cases with onset at start of study
maxIP=12;

% Set maximum incubation period for cases with onset shortly after
% migration in
maxIP_IM=4;

%% MAKE VECTORS OF EVENT TIMES
% Make STATA month origin for converting months in cleaned data
origin=stata_month(startyr,startmo)-1;
% End time
tmax=stata_month(endyr,endmo)-origin;

% Make indicator variables
INTMIG_OUT=(data.INTMIG_OUT==1); % 1st observations for individuals who internally migrated (relocated to a different HH within the study area)
INTMIG_IN=(data.INTMIG_IN==1); % 2nd observations for individuals who internally migrated 
KothrObs=((data.KA>=data.MIG_OUT&INTMIG_OUT)|(data.KA<data.MIG_IN&INTMIG_IN)); % observations for KA cases who internally migrated in which they did not have KA onset
PothrObs=((data.PKDL>=data.MIG_OUT&INTMIG_OUT)|(data.PKDL<data.MIG_IN&INTMIG_IN)); % observations of PKDL cases who internally migrated in which they did not have PKDL onset 
K02_10=(data.KA_1210&(data.KA>origin|(data.KAYR>=startyr&isnan(data.KA)))&~KothrObs); % KA with onset during study (between start and end times)
actvK=((data.KA<=origin&data.KARX>origin)|(data.KAYR==startyr-1&isnan(data.KARX)))&isnan(data.MIG_IN); % (potentially) active KA at start of study
prevK=(data.KAYR<startyr&~actvK); % KA before the start of the study
IpreEXTIM=(K02_10&data.EXTMIG_IN==1&data.KA<=data.MIG_IN&data.KA>origin+maxIP); % KA onset before or at migration in
EXTIMsoonI=(K02_10&data.EXTMIG_IN==1&data.KA>=data.MIG_IN&data.KA-data.MIG_IN<=maxIP_IM&data.KA>origin+maxIP); % KA onset within 6 months of migration in
IpreINTIM=(data.KA<=data.MIG_IN&INTMIG_IN&~prevK); % KA onset before or at internal migration in
PpreINTIM=(data.PKDL<=data.MIG_IN&INTMIG_IN); % PKDL onset before or at internal migration in
PpreEXTIM=(isnan(data.KA)&data.PKDL<=data.MIG_IN&data.EXTMIG_IN==1); % PKDL onset (without prior KA) before or at external migration in
RXF=(data.ALLRXFAIL==1&(data.MOS_RX_NEW_SX==0|(data.PKDL-data.KARX==1))&~KothrObs); % KA cases who suffered treatment failure (include CHMP78102 who had PKDL onset 1 month after KA treatment)
REL=(data.ALLRXFAIL==1&~(data.MOS_RX_NEW_SX==0|(data.PKDL-data.KARX==1))&~KothrObs&~(data.KARX+data.MOS_RX_NEW_SX>data.MIG_OUT)); % KA cases who suffered relapse (exclude PJNA58204 who suffered relapse after migrating out of study area)
RXP=(data.RXD_PKDL==1); % PKDL cases who were treated

% Make event time vectors
tI=data.KA-origin; % KA onset
tR=data.KARX+durRX-origin; % KA recovery
% Add 1 month of infectiousness for treatment failure KA cases
tR(RXF)=tR(RXF)+1;
tRL=data.KARX+data.MOS_RX_NEW_SX-origin; % KA relapse
tRL(RXF)=NaN; % overwrite "relapse" times for treatment failures
tRLR=NaN(n,1); % relapse recovery
tP=data.PKDL-origin; % PKDL onset
tRP=data.PKDL+data.PKDL_DUR-origin; % PKDL recovery
% Overwrite resolution times with treatment times for treated PKDL cases.
tRP(RXP)=max(data.PKRX(RXP),data.PKRX2(RXP))-origin; % use 2nd treatment time for case with 2 PKDL treatments
tB=data.DOB-origin; % birth
tD=data.DEATH-origin; % death
tIM=data.MIG_IN-origin; % immigration
tEM=data.MIG_OUT-origin; % emigration
 
% Make index vectors for different groups of individuals
I=find(K02_10); % KA cases with onset between start month and end month
OR=find(K02_10&~isnan(tI)&(~isnan(tR)|~isnan(tD))); % KA cases with onset and treatment/death time
NONR=find(K02_10&isnan(tI)&isnan(tR)); % KA cases missing onset and treatment time
ONR=find(K02_10&~isnan(tI)&isnan(tR)&isnan(tD)); % KA cases with onset time but no treatment time who didn't die from KA
RNO=find(K02_10&isnan(tI)&~isnan(tR)); % KA cases with treatment time but no onset time
NO=[NONR;RNO]; % KA cases missing onset time
RLO=find(REL&~isnan(tRL)); % KA relapse cases with relapse times
RLNO=find(REL&isnan(tRL)); % KA relapse cases with missing relapse time
RL=[RLO;RLNO];
P=find(~isnan(tP)&~PothrObs); % PKDL cases
IandP=sort(intersect(I,P)); % KA cases who developed PKDL
PIA=setdiff(setdiff(P,I),find(PpreEXTIM)); % PKDL cases w/o KA during study
PI=intersect(PIA,union(find(actvK|prevK),find(tP>tIM&tI<tIM))); % PKDL with KA onset before start of study or before immigration
PA=setdiff(PIA,PI); % PKDL without prior KA (N.B. none of these PKDL cases internally migrated, so don't need to worry about internal migration with asymptomatic infection/PKDL for them)
B=find(tB>0); % exclude people estimated to have been born before START time
D=find(~isnan(tD)); % individuals who died
DpreR=find(~prevK&~isnan(tI)&isnan(tR)&~isnan(tD)); % KA cases during the study who died before treatment
RpreD=setdiff(I,DpreR); % KA cases who were treated before dying
IM=find(~isnan(tIM)); % individuals who migrated into a household (HH) in the study area
EM=find(~isnan(tEM)); % individuals who migrated out of a HH in the study area
IMI=find(KothrObs&tI<tIM&tR>tIM); % active KA cases who internally migrated before being treated
IMP=find(PothrObs&tP<tIM&tRP>tIM|PpreEXTIM); % active PKDL cases who internally migrated before being treated
% Create index vectors for 1st and 2nd observations of internal migrators 
% (complicated definition to allow for subsetting of the data and to ensure
% 1st and 2nd obs are matched in IM_OUT and IM_IN)
IM_OUT=find(ismember(data.RESP_ID,data.ORIG_ID(INTMIG_IN))&INTMIG_OUT); % observation in 1st HH
IM_IN=find(ismember(data.ORIG_ID,data.RESP_ID(INTMIG_OUT))&INTMIG_IN); % observation in 2nd HH
% KA cases who (potentially) initially had active KA
AOR=find(actvK&~isnan(tI)&~isnan(tR)); % cases with active KA
AONR=find(actvK&~isnan(tI)&isnan(tR)); % cases with KA onset before start of study but missing recovery time
ANONR=find(actvK&isnan(tI)&isnan(tR)); % cases with missing KA onset and recovery time but onset in year preceding study start
A=sort([AOR;AONR;ANONR]);
IPNIA=[I;PI;A;IMI;IMP]; % all KA and PKDL cases with symptoms during the study

% Numbers of individuals in different groups
nI=numel(I); % KA cases
nOR=numel(OR); % KA cases with both onset and treatment times
nNONR=numel(NONR); % KA cases w/o onset or treatment times
nONR=numel(ONR); % KA cases w/ onset time but no treatment time
nRNO=numel(RNO); % KA cases w/ treatment time but no onset time
nNO=numel(NO); % KA cases w/ no onset time
nRL=numel(RL); % KA relapse cases
nRLO=numel(RLO); % KA relapse cases w/ relapse time
nRLNO=numel(RLNO); % KA relapse cases w/o relapse time
nIandP=numel(IandP); % KA cases w/ onset between startyr and endyr and later PKDL
nIMI=numel(IMI); % internal immigrants with KA at time of migration
nIMP=numel(IMP); % internal immigrants with PKDL at time of migration
nPI=numel(PI); % PKDL cases w/ KA onset before startyr or before immigration
nPA=numel(PA); % PKDL cases w/o prior KA
nA=numel(A); % KA cases who may have had active KA at start of study
nAONR=numel(AONR); % potentially active KA cases without treatment time
nANONR=numel(ANONR); % potentially active KA cases without onset or treatment time
nIPNIA=numel(IPNIA); % potential infection sources

%% DRAW INITIAL MISSING KA ONSET AND TREATMENT TIMES
% Create vectors of lower and upper bounds for onset month
tIlb=NaN(n,1);
tIlb(I)=max(stata_month(data.KAYR(I),1)-origin,tB(I)+2); % +2 because individuals can only be infected 1 month after birth or immigration
tIub=NaN(n,1);
tIub(I)=min(min(min(stata_month(data.KAYR(I),12)-origin,tR(I)-1),tP(I)-2),tD(I)-1);
% Fit negative binomial distribution to onset-to-treatment (OT) times
[r0,p0]=FitOTdistn(tI,tR); % use all cases with onset and treatment times (not only those with onset in 2002-2010) as this distn is mostly used to impute missing times for active KA cases at start of study

% For individuals with missing onset, diagnosis and treatment times draw
% onset time at random from onset year
for i=1:nNONR
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
    tI(RNO(i))=randi([tIlb(RNO(i)),tIub(RNO(i))],1);
end

% Create vectors of lower and upper bounds for onset month for active KA
% cases at start of study
tIlbA=NaN(n,1);
tIlbA(A)=max(stata_month(data.KAYR(A),1)-origin,tB(A)+2);
tIubA=NaN(n,1);
tIubA(A)=min(min(min(stata_month(data.KAYR(A),12)-origin,tR(A)-1),tP(A)-2),tD(A)-1);

for i=1:nANONR
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

%% DRAW INITIAL KA RELAPSE AND RELAPSE TREATMENT TIMES
p4=mle(tRL(RLO)-tR(RLO)-1,'distribution','geo');
% Relapse times for relapsers with missing relapse time
for i=1:nRLNO
    j=RLNO(i);
    while isnan(tRL(j)) || tRL(j)>min(min(min(tEM(j),tP(j)),tD(j))-1,tmax)-1
        tRL(j)=tR(j)+geornd(p4)+1;
    end
end
% Relapse treatment times for all relapsers
for i=1:nRL
    j=RL(i);
    while isnan(tRLR(j)) || tRLR(j)>min(min(min(tEM(j),tP(j)),tD(j))-1,tmax)
        tRLR(j)=tRL(j)+nbinrnd(r0,p0)+1;
    end
end

%% DRAW INITIAL PRE-SYMPTOMATIC INFECTION TIMES
% Draw infection times from negative binomial distribution with parameters
% r1 and p1
tE=NaN(n,1); % initialise infection time vector
tE(I)=-Inf; % set infection times so that they will all be updated below
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
I1=find(tI>maxIP&~IpreEXTIM&~EXTIMsoonI&~KothrObs);

%% INFECTION PRESSURE FROM KA AND PKDL CASES
% Calculate distances between HHs
[dHH,ia,ib]=CalcHHDists(data);
nHH=numel(ia); % number of HHs
f=histcounts(ib,1:nHH+1)'; % number of individuals in each HH
f=repmat(f,1,length(f)); % a matrix containing column vectors(occupancy count per household) duplicated
ftransp = f';
if strcmp(typ,'Cauchy') % calculate square of distance matrix if using Cauchy kernel
    dHHsqrd=dHH.^2;
else
    dHHsqrd=[];
end
% Calculate spatial kernel 
[KHH,K0old]=Knl_fast(dHH,dHHsqrd,alpha0,beta0,n,nHH,ib,f,ftransp,1:n);

d0=speye(nHH); % sparse identity matrix for additional within-HH transmission
% Calculate HH-level transmission rate matrix
rateHHA=beta0*KHH+delta0*d0;
% Subset and expand to get the transmission rate from KA and PKDL cases to all HH
rateHH=rateHHA(:,ib(IPNIA));

% Infectiousness over time
% Make vector of lesion-specific infectiousnesses for PKDL cases
hv=NaN(n,1);
hv(strcmp(data.TYPE,'MAC_PAP'))=h1; % macular
hv(strcmp(data.TYPE,'PLQ'))=h2; % plaque
hv(strcmp(data.TYPE,'NOD'))=h3; % nodular
hv(~isnan(tP)&strcmp(data.TYPE,''))=hmssng; % unexamined

% Construct infectiousness matrix (rows = individuals, columns = times)
h=zeros(nIPNIA,tmax); % initialise infectiousness matrix
% PKDL cases w/ prior KA
for i=1:nIandP
    j=IandP(i);
    h(I==j,max(tP(j),tIM(j))+1:min(min(tRP(j),tEM(j)),tmax))=hv(j);
end
% KA cases
for i=1:nI
    j=I(i);
    h(i,max(max(0,tIM(j)),tE(j))+1:tI(j))=h0; % infectiousness from month of infection up to month before symptom onset
    h(i,max(tIM(j),tI(j))+1:min(tEM(j),tRorD(j)))=1; % infectiousness up to month before month of treatment or death
end
% PKDL w/ KA onset before start of study
for i=1:nPI
    j=PI(i);
    h(nI+i,tP(j)+1:min(min(tRP(j),tEM(j)),tmax))=hv(j);
end
% KA relapse cases - assume individuals return to full infectiousness upon relapse until re-treatment
for i=1:nRL
    j=RL(i);
    h(IPNIA==j,tRL(j)+1:min(min(tRLR(j),tEM(j)),tmax))=1;
end
% initially active KA cases
for i=1:nA
    j=A(i);
    h(nI+nPI+i,1:tRorD(j))=1;
end
% KA cases who internally migrated while they had symptoms
for i=1:nIMI
    j=IMI(i);
    h(nI+nPI+nA+i,tIM(j)+1:tRorD(j))=1;
end
% PKDL cases who internally migrated while they had symptoms
for i=1:nIMP
    j=IMP(i);
    h(nI+nPI+nA+nIMI+i,tIM(j)+1:min(min(tRP(j),tEM(j)),tmax))=hv(j);
end
% Convert infectiousness matrix to sparse matrix
h=sparse(h);

% Calculate HH-level infection pressure from KA and PKDL cases
lambdaHHI=rateHH*h; % HH-level infection pressure
lambdaI=lambdaHHI(ib,:); % expand to individual-level infection pressure

%% DRAW INITIAL ASYMPTOMATIC INFECTION AND RECOVERY TIMES
% Index vector for individuals who were either susceptible or
% asymptomatically infected by the end of the study
SusA=setdiff((1:n)',union(union(IPNIA,PA),find(actvK|prevK|IpreEXTIM|EXTIMsoonI|IpreINTIM|PpreINTIM|PpreEXTIM|KothrObs|PothrObs)));
nSusA=numel(SusA);
% Find age of all individuals (in months) at start of study (t=0). Set to 0
% if individual was not yet born
age=max(-tB,0);
% Get range of times for which individual was present in study area 
rng=[max(max(tB,tIM)+1,0) min(min(tEM,tD),tmax+1)];
% Index vector for individuals present at (born or imigrated before) start of study
Pres0=find(rng(:,1)==0);
% figure; histogram(rng(:,1),'BinM','int'); hold on; histogram(rng(:,2),'BinM','int'); hold off
% Index vector for individuals who were either initially susceptible or
% asymptomatically infected
SusA0=intersect(SusA,Pres0);

% Calculate probabilities of each individual being initially susceptible or 
% actively asymptomatically infected according to initial status model
prob0=ProbInitStatus(age,lambda0,p2);
% Append probabilities of individuals being initially recovered from 
% previous asymptomatic infection
prob0=[prob0,1-sum(prob0,2)];
% Overwrite probability of being initially susceptible for individuals who
% were not born or had not yet immigrated as 1, and other probabilities as 0
prob0(rng(:,1)>0,1)=1;
prob0(rng(:,1)>0,2:3)=0;

% Draw initial statuses for non-symptomatic individuals using calculated
% probabilities
Stat0=NaN(n,1);
for i=1:numel(SusA0)
    j=SusA0(i);
    Stat0(j)=randsample(1:3,1,true,prob0(j,:));
end
% Make index vectors for initial infection status
actvA=find(Stat0==2); % actively asymptomatically infected
nactvA=numel(actvA);
prevA=find(Stat0==3); % previously asymptomatically infected

% Make index vector for 2nd observations for internal migrators who were
% intially previously/actively asymptomatically infected
IM_INprevAactvA=IM_IN(ismember(IM_OUT,[prevA;actvA]));
% Remove individuals initially recovered from asymptomatic infection or 
% actively asymptomatically infected and their 2nd observations (if they 
% internally migrated) from index of susceptible individuals  
Sus=setdiff(SusA,[union(prevA,actvA);IM_INprevAactvA]);

% For internal migrators, randomly pick one observation to treat as
% susceptible (as they can only *become* infected during one observation)
IM_OUT_IN=[IM_OUT,IM_IN]; % paired indices for 1st and 2nd observations for internal migrators
IM_OUTSus=intersect(IM_OUT,Sus); % 1st observations for susceptible internal migrators 
nIM_OUTSus=numel(IM_OUTSus); % number of susceptible internal migrators
rmIM=NaN(nIM_OUTSus,1); % observation to remove
kpIM=NaN(nIM_OUTSus,1); % observation to keep
for i=1:nIM_OUTSus
    rmIM(i)=IM_OUT_IN(IM_OUT==IM_OUTSus(i),randi(2));
    kpIM(i)=IM_OUT_IN(IM_OUT==IM_OUTSus(i),IM_OUT_IN(IM_OUT==IM_OUTSus(i),:)~=rmIM(i));
end
Sus=setdiff(Sus,rmIM); % remove selected observations from index of susceptibles
nSus=numel(Sus);
% Make indicator matrix of times individuals were present in study area 
% (rows = individuals, columns = times)
rngm=ones(n,tmax);
for i=1:n
    if rng(i,1)>1 || rng(i,2)<tmax+1 % if individual entered after t=1 or left before t=tmax+1
        rngm(i,[1:rng(i,1)-1,rng(i,2):tmax])=0; % zero out entries in matrix corresponding to when they weren't present
    end
end
% Calculate cumulative infection pressure on individuals from KA and PKDL
% cases up to each time point
cum_lambdaI=cumsum(lambdaI.*rngm,2);
cum_lambdaI=[zeros(n,1),cum_lambdaI(:,1:end-1)];
% Calculate the probability of individuals being asymptomatically infected
% at each time point, or avoiding asymptomatic infection for the duration 
% of the study, based on the infection pressure from KA and PKDL cases
probA1=[exp(-cum_lambdaI).*(1-exp(-(1-pI0)*lambdaI.*rngm)),exp(-sum(lambdaI.*rngm,2))];
probA1=bsxfun(@rdivide,probA1,sum(probA1,2)); % normalise the probabilities
probA=[1-prob0(:,1),bsxfun(@times,prob0(:,1),probA1)]; % multiply by the probability of being initially susceptible and prepend probability of initially having been asymptomatically infected
% Calculate probability of having been asymptomatically infected before the
% end of the study
cum_probA=sum(probA(:,1:end-1),2);

% Set initial asymptomatic infection and recovery times
tA=NaN(n,1);
tRA=NaN(n,1);
% Set infection and recovery times for individuals initially previously
% asymptomatically infected to 0
tA(prevA)=0;
tRA(prevA)=0;
% Set infection time for individuals initially actively asymptomatically
% infected to 0
tA(actvA)=0;
% Draw recovery times for individuals initially actively asymptomatically
% infected
for i=1:numel(actvA)
    j=actvA(i); % get index of ith initially actively asymptomatically infected individual
    if ismember(j,IM_OUT) % if an internal migrator
        j1=IM_IN(IM_OUT==j); % get index of 2nd observation
        probAIP=[geopdf(0:rng(j1,2)-2,p2),1-geocdf(rng(j1,2)-2,p2)]; % calculate probabilities of possible asymptomatic infection durations
        fwdRA=[1:rng(j1,2)-1,tmax+1]; % make vector of possible asymptomatic recovery times
    else % if not an internal migrator
        probAIP=[geopdf(0:rng(j,2)-2,p2),1-geocdf(rng(j,2)-2,p2)]; % calculate probabilities of possible asymptomatic infection durations
        fwdRA=[1:rng(j,2)-1,tmax+1]; % make vector of possible asymptomatic recovery times
    end    
    tRA(j)=fwdRA(randsample(numel(fwdRA),1,true,probAIP)); % draw recovery time
end
% Set infection and recovery times for 2nd observations for internal 
% migrators who were initially previously or actively asymptomatically 
% infected to dummy time of tmax+2
tA(IM_INprevAactvA)=tmax+2;
tRA(IM_INprevAactvA)=tmax+2;
% Make index vectors for 1st and 2nd observations of internal migrators 
% with active asymptomatic infection at start of study
IM_OUTactvA=IM_OUT(ismember(IM_OUT,actvA));
IM_INactvA=IM_IN(ismember(IM_OUT,actvA));
% Make index vectors for internal migrators with active asymptomatic 
% infection at the start of the study who recovered before internal
% migration (during their 1st obs) and those who recovered after internal 
% migration (during their 2nd obs)
RAobs1actvA=IM_OUTactvA(tRA(IM_OUTactvA)>rng(IM_OUTactvA,2)-1);
RAobs2actvA=IM_INactvA(tRA(IM_OUTactvA)>rng(IM_OUTactvA,2)-1);

% Calculate the number of asymptomatic infections without subsequent PKDL
% during the study using the (initial) proportion of infections that lead
% to KA (pI0), the observed number of KA cases (nI) and the observed number
% of PKDL cases w/o prior KA
nAsx=min(round((1-pI0)/pI0*nI)-nPA,nSus); % limit the number of asymptomatic infections to the size of the susceptible population regardless of the value of pI0
% Draw the individuals who get asymptomatically infected according to the
% probability of having been asymptomatically infected before the end of
% the study (i.e. using a weighted draw)
Asx=Sus(datasample(1:nSus,nAsx,'Replace',false,'Weights',cum_probA(Sus)));
% Make index vector of individuals still susceptible at end of study
Susend=[setdiff(Sus,Asx);rmIM(~ismember(kpIM,Asx));IM_OUT(ismember(IM_IN,Asx))];
tA(Susend)=tmax+1;
tRA(Susend)=tmax+1;
% Set asymptomatic infection and recovery time of 2nd observation of
% internal migrators asymptomatically infected during 1st observation to
% dummy time of tmax+2
tA(IM_IN(ismember(IM_OUT,Asx)))=tmax+2;
tRA(IM_IN(ismember(IM_OUT,Asx)))=tmax+2;
% Draw asymptomatic infection and recovery times according to probability
% of asymptomatic infection from KA and PKDL cases
for i=1:nAsx
    j=Asx(i);
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
% Make index vectors for 1st and 2nd observations of internal migrators who
% get asymptomatically infected
IM_OUTAsx=IM_OUT(ismember(IM_OUT,Asx));
IM_INAsx=IM_IN(ismember(IM_OUT,Asx));
% Make index vectors for internal migrators asymptomatically infected 
% during the study who recovered before internal migration (during their 
% 1st obs) and those who recovered after internal migration (during their 
% 2nd obs)
RAobs1=IM_OUTAsx(tRA(IM_OUTAsx)>rng(IM_OUTAsx,2)-1);
RAobs2=IM_INAsx(tRA(IM_OUTAsx)>rng(IM_OUTAsx,2)-1);

% % Fit negative binomial distribution to dormant infection period
% % distribution
% RP=tP-tR;
% pars1=nbinfit(RP(RP>=0));
% r3=pars1(1);
% p3=pars1(2);
% Load parameters for negative binomial distribution fitted to dormant
% infection period for full dataset
load('FixedParams','r3'); 
load('FixedParams','p3'); 
% Use this to draw asymptomatic recovery times and infection times for PKDL
% cases w/o prior KA
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

% Make index vector of individuals actively asymptomatically infected at
% the start of the study or those asymptomatically infected during the
% study
A1=[Asx;actvA];
nA1=numel(A1);

%% INFECTIOUSNESS
% Make infectiousness matrix for asymptomatics
hA=zeros(n,tmax);
for i=1:nA1
    j=A1(i);
    hA(j,tA(j)+1:min(min(min(tRA(j),tEM(j)),tD(j)),tmax))=h40;
end
% Add infectiousness to 2nd observations for internal migrators with active
% asymptomatic infection at the start of the study who recovered during
% their 2nd observation
for i=1:numel(RAobs2actvA)
    j=RAobs1actvA(i);
    j1=RAobs2actvA(i);
    hA(j1,rng(j,2)+1:min(min(tRA(j),rng(j1,2)),tmax))=h40;
end 
% Add infectiousness to 2nd observations for internal migrators 
% asymptomatically infected during the study who recovered during their 2nd 
% observation
for i=1:numel(RAobs2)
    j=RAobs1(i);
    j1=RAobs2(i);
    hA(j1,rng(j,2)+1:min(min(tRA(j),rng(j1,2)),tmax))=h40;
end
% Convert to a sparse matrix
hA=sparse(hA);

% Make infectiousness matrix for PKDL cases w/o prior KA 
hPA=zeros(nPA,tmax);
for i=1:nPA
    j=PA(i);
    hPA(i,max(1,tA(j)+1):tRA(j))=h40;
    hPA(i,tP(j)+1:min(min(tRP(j),tEM(j)),tmax))=hv(j);
end
% figure; spy(h); axis square
% figure; spy(hA); axis square
% figure; spy(hPA); axis square

%% TOTAL INFECTION PRESSURE
% Calculate HH-level infection pressure from all asymptomatic individuals
% w/o subsequent PKDL
lambdaHHA=rateHHA(:,ib([A1;RAobs2actvA;RAobs2]))*hA([A1;RAobs2actvA;RAobs2],:);
% Calculate HH-level infection pressure from PKDL cases w/o prior KA
rateHHPA=rateHHA(:,ib(PA));
lambdaHHPA=rateHHPA*hPA;
% Calculate total HH-level infection pressure from all infectious
% individuals
lambdaHH=rateHH*h+lambdaHHPA+lambdaHHA+epsilon0;
lambda=lambdaHH(ib,:); % expand to individual-level infection pressure

%% MAKE STATUS MATRICES
% Make logical matrices for pre-symptomatic infection, asymptomatic 
% infection, birth, death, immigration and emigration times for
% constructing susceptibility matrix
% Pre-symptomatic infection
tEm=false(n,tmax);
tEm((tE(I1)-1)*n+I1)=1;
tEm=sparse(tEm);
% Asymptomatic infection
tAm=false(n,tmax);
% Index vector of individuals asymptomatically infected during study
A2=find(tA>0 & tA<tmax+1);
tAm((tA(A2)-1)*n+A2)=1;
tAm=sparse(tAm);
% figure; spy(tEm); hold on; spy(tAm,'r'); axis square
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

% Make susceptibility status matrix
S=1-max(preB,preIM)-max(max(max(cumsum(tEm,2),cumsum(tAm,2)),cumsum(tDm,2)),cumsum(tEMm,2));
% Remove susceptibility for various groups who were not susceptible when 
% they entered the study or had onset within the initial incubation period 
% window
S(prevK,:)=0; % previous KA cases
S(prevA,:)=0; % previously asymptomatically infected individuals
S(actvA,:)=0; % initially actively asymptomatically infected individuals
S(tI<=maxIP,:)=0; % KA cases with onset before maxIP (also excludes cases with active KA at start of study and a couple of cases with onset in 2002 before immigration)
S(IpreEXTIM,:)=0; % KA cases with onset before or at migration in
S(EXTIMsoonI,:)=0; % KA cases with onset within 6 months of migration in
S(IpreINTIM,:)=0; % KA cases with onset before internal migration in
S(PpreINTIM,:)=0; % PKDL cases with onset before internal migration in
S(PpreEXTIM,:)=0; % PKDL cases (without prior KA) with onset before external migration in
S(prevAPA,:)=0; % previously asymptomatically infected individuals who later developed PKDL
S(actvAPA,:)=0; % initially actively asymptomatically infected individuals who later developed PKDL
S(IM_INprevAactvA,:)=0; % 2nd observations of internal migrators asymptomatically infected before the start of the study
S(IM_IN(ismember(IM_OUT,Asx)),:)=0; % 2nd observations of internal migrators asymptomatically infected during 1st observation

% figure; histogram(S)

% % Checks
% figure; spy(S(setdiff(I,IandP),:)); axis square; pause(2); hold on
% spy(h(~ismember(I,IandP),:),'r'); axis square; hold off
% figure; spy(S(IandP,:)); axis square; pause(2); hold on
% spy(h(ismember(I,IandP),:),'r'); axis square; hold off
% figure; spy(S([IMI-1,IMI,IMP(1)-1,IMP(1),IMP(2)-1,IMP(2)],:)); axis square; pause(2); hold on
% spy(h([find(I==IMI-1),end-2,find(I==IMP(1)-1),end-1,nI+find(PNI==IMP(2)-1),end],:),'r'); axis square; hold off

% Make index vector for individuals present and susceptible at start of
% study
S0=intersect(union(Sus,find(S(:,1))),Pres0);

%% SET UP MCMC
% Set MCMC parameters
burnin=round(niters/10); % burn-in
nu=numel(u); % number of parameters to update
pname={'beta','alpha','epsilon','delta','lambda_0','h_4','p_I'}; % parameter names
np=numel(pname); % number of parameters that can be updated

% Set prior distributions for parameters
prior_shape=ones(1,6);
prior_scale=[10,100,1,1,0.01,h40];
priorpdf={'gam','gam','gam','gam','gam','gam','beta'};
priorp=cell(1,np);
for i=1:6
    priorp{i}=[prior_shape(i),prior_scale(i)];
end
priorp{7}=[15,85]; %[1,1];

% Set initial proposal covariance matrix for transmission parameter block update
ppvar0=diag([0.01,100,1e-7,4e-4,4e-8,1e-4,1e-4]/25);

% Set initial value for running estimate of covariance matrix
ppvar=ppvar0;

% Make vector for storing scale factor for proposal covariance matrix for block update
c=NaN(niters+1,1);
c(1)=1; % initial scale factor

% Set number of pre-symptomatic and asymptomatic infection times to propose new values for per iteration
nEmoves=round(nOR/5); % number of pre-symptomatic infection times
% Randomly select which KA cases to update infection times for in each iteration
pick=zeros(nEmoves+nNONR+nONR+nRNO,niters);
for i=1:niters 
pick(:,i)=[OR(randperm(nOR,nEmoves));NONR;ONR;RNO]; % N.B. randperm rather than randi to avoid possibility of changing same infection time twice in one step
end
nAmoves=round(nAsx/5);%round(nAsx/10);%round((round((1-pI)/pI*nI)-nPA)/5);%round(nAsx/3);%round(nAsx/10);%round(nAsx/20);%
pickA=zeros(nAmoves,niters,'uint16');
% for i=1:niters
% pickA(:,i)=SusA(randperm(nSusA,nAupdts));
% end
% Set maximum allowable asymptomatic infection time move
M=tmax;
% Proposal variance for infection period moves for KA cases missing onset
% and treatment times
ERvar=4; % months

% Initialisation
pold=[beta0,alpha0,epsilon0,delta0,lambda0,h40,pI0]; % initial transmission parameter values
ppmean=pold; % initial proposal distribution means
plb=[zeros(1,6),0]; % lower bounds for parameters
pub=[Inf(1,6),1]; % upper bounds for parameters
p1new=p10; % initial incubation period distribution parameter value
% Initial incubation periods
IPold=NaN(n,1);
IPold(I)=tI(I)-tE(I);
% Initial values of log-likelihood terms
LL1old=L1(S,lambda);
LL2old=L2(lambda,pI0,tEm);
LL3old=L3(IPold(I),r1,p10);
LL4old=L4(lambda,pI0,tAm);
LL5old=L5(age,S0,actvA,prevA,lambda0,pI0,p2);
LL6old=L5(age,S0PA,actvAPA,prevAPA,lambda0,pI0,p2);
 
% Matrices and vectors for saving parameters, log-likelihood values and 
% missing event times
p=zeros(niters+1,np); % transmission parameter values
p(1,:)=pold; % initial transmission parameter values
p1=zeros(niters,1); % p1 values
K0=zeros(niters,1); % spatial kernel normalisation constants
LL=zeros(niters,1); % log-likelihood values
% terms=zeros(niters,5); % individual log-likelihood terms
terms=zeros(niters,6); % individual log-likelihood terms
tEs=zeros(nI,niters,'int8'); % presymptomatic infection times (integer matrix to save memory)
IPs=zeros(nI,niters); % incubation periods
% N.B. Matrices for saving asymptomatic infection and recovery times need 
% to be integer matrices otherwise they use too much memory to be saved
% without thinning
tAs=zeros(n,niters,'uint8'); % asymptomatic infection times
tAs(setdiff((1:n)',[SusA;PA]),:)=tmax+2; % use dummy asymptomatic infection time for symptomatic individuals
tRAs=zeros(n,niters,'uint8'); % asymptomatic recovery times
tRAs(setdiff((1:n)',[SusA;PA]),:)=tmax+2; % use dummy asymptomatic recovery time for symptomatic individuals
tIsNONR=zeros(nNONR,niters); % onset times of KA cases missing onset and treatment times
tRsNONR=zeros(nNONR,niters); % treatment times of KA cases missing onset and treatment times
tIsRNO=zeros(nRNO,niters); % onset times of KA cases missing only onset time
tRsONR=zeros(nONR,niters); % treatment times of KA cases missing only treatment time
tIsANONR=zeros(nANONR,niters); % onset times of potentially initially active KA cases missing onset and treatment times
tRsANONR=zeros(nANONR,niters); % treatment onset times of potentially initially active KA cases missing onset and treatment times
tRsAONR=zeros(nAONR,niters); % recovery times of potentially initially active KA cases missing only treatment time
tRLRsRLO=zeros(nRLO,niters); % relapse treatment times for KA relapse cases with observed relapse time
tRLsRLNO=zeros(nRLNO,niters); % relapse times for KA relapse cases with missing relapse times
tRLRsRLNO=zeros(nRLNO,niters); % relase treatment times for KA relapse cases with missing relapse times
 
% Initialise new parameter values, transmission rate matrices, event time vectors and infectiousness matrices
pnew=pold;
rateHH_new=rateHH;
rateHHA_new=rateHHA;
rateHHPA_new=rateHHPA;
lambdaHH_new=lambdaHH;
lambdaHHA_new=lambdaHHA;
lambdaHHPA_new=lambdaHHPA;
lambda_new=lambda;
lambda_mean=lambda; % running estimate of individual-level infection pressure (mean over the chain)
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
A1new=A1;
actvAnew=actvA;
prevAnew=prevA;
S0new=S0;
actvAPAnew=actvAPA;
prevAPAnew=prevAPA;
S0PAnew=S0PA;
RAobs2actvAnew=RAobs2actvA;
RAobs2new=RAobs2;

% Initialise acceptance and rejection counts and rates
% Transmission parameter updates
acc_p=0;
rej_p=0;
acc_rate_p=0;
% Pre-symptomatic infection time updates
acc_E=0;
rej_E=0;
acc_rate_E=0;
% Asymptomatic infection time moves from t=T+1 to t<T+1
acc_Aadd=0;
rej_Aadd=0;
acc_rate_Aadd=0;
% Asymptomatic infection time moves from t=0 to t>0
acc_Arem=0;
rej_Arem=0;
acc_rate_Arem=0;
% Asymptomatic infection time moves from t in [1,T] to t' in [1,T]
acc_Amov=0;
rej_Amov=0;
acc_rate_Amov=0;
% Individual-level asymptomatic infection time moves
acc_add=zeros(n,1);
rej_add=zeros(n,1);
acc_rem=zeros(n,1);
rej_rem=zeros(n,1);
acc_mov=zeros(n,1);
rej_mov=zeros(n,1);
% KA cases with missing onset times
acc_I=0;
rej_I=0;
acc_rate_I=0;
% KA cases missing both onset and treatment times 
acc_ERmove=0;
rej_ERmove=0;
acc_rate_ERmove=0;
% KA cases missing treatment time
acc_R=0;
rej_R=0;
acc_rate_R=0;
% Potentially initially active KA cases with missing onset and treatment
% times
acc_AIRmove=0;
rej_AIRmove=0;
acc_rate_AIRmove=0;
% Potentially initially active KA cases with missing treatment time
acc_AR=0;
rej_AR=0;
acc_rate_AR=0;
% PKDL cases w/o prior KA
acc_PA=0;
rej_PA=0;
acc_rate_PA=0;
% Relapse cases with missing relapse and relapse treatment time
acc_RLNO=0;
rej_RLNO=0;
acc_rate_RLNO=0;
% Relapse cases with missing treatment time
acc_RLO=0;
rej_RLO=0;
acc_rate_RLO=0;

% Parameters for plotting output
nbins=50;
scrnsz=get(0,'ScreenSize');

%% MCMC LOOP
for k=1:niters
    %% UPDATE TRANSMISSION PARAMETERS USING ADAPTIVE RANDOM WALK METROPOLIS-HASTINGS    
    pnew(u)=mvnrnd(pold(u),c(k)^2*2.38^2*ppvar(u,u)/nu);
    
    if all(pnew(u)>plb(u) & pnew(u)<pub(u)) % check if prior probability is non-zero before calculating log-likelihood
        % Calculate ratio of prior probabilities for new and old value
        q=0;
        for i=u % each parameter that is to be updated
            if numel(priorp{i})==1 % number of hyperparameters in prior distribution is 1
                q=q+log(pdf(priorpdf{i},pnew(i),priorp{i}))-log(pdf(priorpdf{i},pold(i),priorp{i}));
            elseif numel(priorp{i})==2 % number of hyperparameters in prior distribution is 2
                q=q+log(pdf(priorpdf{i},pnew(i),priorp{i}(1),priorp{i}(2)))-log(pdf(priorpdf{i},pold(i),priorp{i}(1),priorp{i}(2)));
            end
        end
        
        % Calculate new infection pressure
        [KHH_new,K0new]=Knl_fast(dHH,dHHsqrd,pnew(2),pnew(1),n,nHH,ib,f,ftransp,1:n); % update spatial kernel
        rateHHA_new=pnew(1)*KHH_new; % new HH-level transmission rate
        if pnew(4)~=0 % additional within-HH transmission rate is non-zero
            rateHHA_new=rateHHA_new+pnew(4)*d0; % add additional within-HH contribution
        end
        rateHH_new=rateHHA_new(:,ib(IPNIA)); % new transmission rate from each KA and PKDL case to each HH
        rateHHPA_new=rateHHA_new(:,ib(PA)); % new transmission rate from each PKDL case w/o prior KA to each HH
        if ismember(6,u) % relative asymptomatic infectiousness is being estimated
            hAnew=pnew(6)/pold(6)*hA; % rescale asymptomatic infectiousness matrix
        end
        lambdaHHA_new=rateHHA_new(:,ib([A1;RAobs2actvA;RAobs2]))*hAnew([A1;RAobs2actvA;RAobs2],:); % new infection pressure on HHs from all asymptomatic individuals
        lambdaHHPA_new=rateHHPA_new*hPA; % new infection pressure on HHs from PKDL cases w/o prior KA
        lambdaHH_new=rateHH_new*h+lambdaHHPA_new+lambdaHHA_new+pnew(3); % new total infection pressure on HHs
        lambda_new=lambdaHH_new(ib,:); % new individual-level infection pressure

        % Calculate new log-likelihood terms
        LL1new=L1(S,lambda_new);
        LL2new=L2(lambda_new,pnew(7),tEm);
        LL4new=L4(lambda_new,pnew(7),tAm);
        LL5new=L5(age,S0,actvA,prevA,pnew(5),pnew(7),p2);
        LL6new=L5(age,S0PA,actvAPA,prevAPA,pnew(5),pnew(7),p2);
        LLnew=LL1new+LL2new+LL4new+LL5new+LL6new;
        LLold=LL1old+LL2old+LL4old+LL5old+LL6old;

        % Calculate Metropolis-Hastings acceptance probability
        log_ap=LLnew-LLold+q;

        if log_ap > log(rand) % accept new values if acceptance probability > rand
            pold=pnew; % keep new parameter values
            lambda=lambda_new; lambdaHH=lambdaHH_new; rateHH=rateHH_new; 
            rateHHA=rateHHA_new; rateHHPA=rateHHPA_new; KHH=KHH_new; K0old=K0new;
            lambdaHHA=lambdaHHA_new; lambdaHHPA=lambdaHHPA_new; hA=hAnew;
            LL1old=LL1new; LL2old=LL2new;
            LL4old=LL4new; LL5old=LL5new; LL6old=LL6new;
            acc_p=acc_p+1; % add to acceptance
            c(k+1)=(1+100/(100+k))*c(k); % increase proposal variance scaling constant
        else % else reject and keep old values
            lambda_new=lambda; lambdaHH_new=lambdaHH;
            lambdaHHA_new=lambdaHHA; lambdaHHPA_new=lambdaHHPA; hAnew=hA;
            rej_p=rej_p+1; % add to rejection
            c(k+1)=(1+100/(100+k))^(0.234/(0.234-1))*c(k); % shrink proposal variance scaling
        end
    else % if prior probability is 0, reject immediately
        rej_p=rej_p+1;
        c(k+1)=(1+100/(100+k))^(0.234/(0.234-1))*c(k); % shrink proposal variance scaling
    end
    
    %% GIBBS UPDATE SUCCESS PROBABILITY PARAMETER FOR INCUBATION PERIOD DISTRIBUTION
    p1new=betarnd(a+r1*nI,b+sum(IPold(I))-nI);
 
    %% UPDATE PRE-SYMPTOMATIC INFECTION TIMES    
    for i=1:size(pick,1)
        j=pick(i,k); % get index of infection time to update
        tEnew(j)=tI(j)-(nbinrnd(r1,p1new)+1); % draw new infection time by drawing new incubation period and subtracting it from onset time
        
        if tEnew(j)>=tB(j)+1 && ~(tI(j)>maxIP && tEnew(j)<1) && ~(tI(j)-tIM(j)-1>maxIP_IM && tEnew(j)<tIM(j)+1) % calculate log-likelihood if new infection time is after birth, not before start of study if onset is after max_IP, and not before immigration if onset is more than maxIP_IM months after immigration
            tEj=tE(j);
            tEjnew=tEnew(j);
            tIMj=tIM(j);
            IPnew(j)=tI(j)-tEjnew; % new incubation period
            q=log(nbinpdf(IPold(j)-1,r1,p1new))-log(nbinpdf(IPnew(j)-1,r1,p1new)); % proposal ratio for infection time update
            
            % Update infection time and susceptibility matrices
            if tI(j)>maxIP && ~IpreEXTIM(j) && ~EXTIMsoonI(j) % only for cases with onset after initial window and after entry
                tEmnew(j,tEj)=0; % remove old infection time
                tEmnew(j,tEjnew)=1; % add new infection time
                Snew(j,tEjnew:tEj-1)=0; % remove old susceptible times if new infection time is earlier (from month indvdl is infctd up to but not incl. month before old infection time)
                Snew(j,tEj:tEjnew-1)=1; % add new susceptible times if new infection time is later (up to but not incl. month indvdl is infctd)
            end
            
            % Update infectiousness matrix
            m=(IPNIA==j); % get index of individual in infectiousness matrix
            entry=max(0,tIMj); % get entry time of individual
            hnew(m,max(entry,tEj)+1:tEjnew)=0; % remove old infectiousness if new infection time is later
            hnew(m,max(entry,tEjnew)+1:tEj)=h0; % add infectiousness if new infection time is earlier
            erlrE=max(entry,min(tEj,tEjnew)); % index of column for earlier infctn time between old and new infection time
            ltrE=max(entry,max(tEj,tEjnew)); % index of column for later infctn time between old and new infection time
            
            % Calculate new infection pressure
            idx=erlrE+1:ltrE; % indices of columns in infection pressure matrix to change
            lambdaHH_new(:,idx)=rateHH*hnew(:,idx)+lambdaHHPA(:,idx)+lambdaHHA(:,idx)+pold(3); % new infection pressure on each HH
            lambda_new(:,idx)=lambdaHH_new(ib,idx); % new infection pressure on each individual

            % Calculate new log-likelihood terms
            LL1new=L1(Snew,lambda_new);
            LL2new=L2(lambda_new,pold(7),tEmnew);
            LL3new=L3(IPnew(I),r1,p1new);
            LL4new=L4(lambda_new,pold(7),tAm);
            LLnew=LL1new+LL2new+LL3new+LL4new;
            LLold=LL1old+LL2old+LL3old+LL4old;
            
            % Calculate Metropolis-Hastings acceptance probability
            log_ap=LLnew-LLold+q;
            
            if log_ap > log(rand) % accept new values if acceptance probability > rand
                IPold(j)=IPnew(j); tE(j)=tEnew(j); S(j,:)=Snew(j,:); 
                tEm(j,:)=tEmnew(j,:); h(m,:)=hnew(m,:); 
                lambda(:,idx)=lambda_new(:,idx); lambdaHH(:,idx)=lambdaHH_new(:,idx);
                LL1old=LL1new; LL2old=LL2new; LL3old=LL3new; 
                LL4old=LL4new; 
                acc_E=acc_E+1;
            else % else reject and keep old values
                IPnew(j)=IPold(j); tEnew(j)=tE(j); Snew(j,:)=S(j,:); 
                tEmnew(j,:)=tEm(j,:); hnew(m,:)=h(m,:); 
                lambda_new(:,idx)=lambda(:,idx); lambdaHH_new(:,idx)=lambdaHH(:,idx); % keep old values, don't change log-likelihood
                rej_E=rej_E+1;
            end
        else % otherwise reject immediately
            tEnew(j)=tE(j);
            rej_E=rej_E+1;
        end
    end
     
    %% UPDATE ASYMPTOMATIC INFECTION TIMES    
    % Recalculate cumulative infection pressure
    cum_lambda=cumsum(lambda_mean.*rngm,2);
    cum_lambda=[zeros(n,1),cum_lambda(:,1:end-1)];
    % Recalculate probabilities of individuals being asymptomatically
    % infected at each time point during the study or avoiding infection
    % for the whole of the study
    probA1=[exp(-cum_lambda).*(1-exp(-(1-pold(7))*lambda_mean.*rngm)),exp(-sum(lambda_mean.*rngm,2))];
    probA1=bsxfun(@rdivide,probA1,sum(probA1,2));
    % Update probabilties of initially being susceptible or actively
    % asypmtomatically infected
    prob0=ProbInitStatus(age,pold(5),p2);
    % Append probabilities of initially being recovered from asymptomatic
    % infection
    prob0=[prob0,1-sum(prob0,2)];
    % Multiply asymptomatic infection probabilities by probability of
    % initially being susceptible and prepend probability of initially
    % having been asymptomatically infected
    probA=[1-prob0(:,1),bsxfun(@times,prob0(:,1),probA1)];
    % Calculate cumulative probability of asymptomatic infection before the
    % end of the study
    cum_probA=sum(probA(:,1:end-1),2);
    
    for i=1:nAmoves
       % Pick non-symptomatic individual with probability proportional to
       % cumulative probability of asymptomatic infection before the end of
       % the study
       j=SusA(sum(rand>=cumsum(cum_probA(SusA))/sum(cum_probA(SusA)))+1);
       mig_out=ismember(j,IM_OUT); % flag for whether j is internal migrator and it's their 1st observation
       mig_in=ismember(j,IM_IN); % flag for whether j is internal migrator and it's their 2nd observation
       pickA(i,k)=j; % save which individual was picked
       
       % Only try to update asymptomatic infection and recovery time of j
       % if they are not an internal migrator or are an internal migrator 
       % who is not already infected under their other observation
       if ~(mig_out || mig_in) || (mig_out && ~ismember(IM_IN(IM_OUT==j),[prevA;A1])) || (mig_in && ~ismember(IM_OUT(IM_IN==j),[prevA;A1]))
       t=tA(j); % get current asymptomatic infection time
       s=tRA(j); % get current asymptomatic recovery time
       idx1=t+1:min(min(s,rng(j,2)),tmax); % index of columns in infectiousness matrix to remove infectiousness from
       if mig_out % if j is internal migrator and this is their 1st observation
           j1=IM_IN(IM_OUT==j); % get index of 2nd observation
           if t<=rng(j,2)-1 && s>rng(j,2)-1 % check if j was infected before internal migration and remained infected after moving
               idx3=rng(j,2)+1:min(min(s,rng(j1,2)),tmax); % index of columns to change in infectiousness matrix for infection during 2nd observation
           end
       end
       if t==0 % currently asymptomatically infected before start of study (N.B. means must be j's 1st observation)
           fwd=[0:rng(j,2)-1,tmax+1]; % times that j's infection time can move to
           tp=fwd(sum(rand>=cumsum(probA(j,fwd+1))/sum(probA(j,fwd+1)))+1); % draw new infection time from fwd using current estimated probability of asymptomatic infection at each of those times
           if tp==0 % if new infection time is before start of study
               bck=[0:rng(j,2)-1,tmax+1]; % possible infection times for reverse move
               if rand<prob0(j,3)/sum(prob0(j,2:3)) % recovered before start of study 
                   sp=0; % new recovery time is 0
                   if s==0 % if j is not already in asymptomatic infectiousness matrix as they recovered by start of study
                       idx=[]; % no change required
                       q=0; % proposal ratio is 0
                   else % if j is already in asymptomatic infectiousness matrix as they recovered after start of study
                       prevAnew=[prevA;j]; % add them to index of initially previously asymptomatically infected individuals
                       actvAnew=actvA(actvA~=j); % remove them from index of initially active asymptomatics
                       A1new=A1(A1~=j); % remove them from index of initially active asymptomatics and individuals asymptomatically infected during the study
                       idx=idx1;
                       hAnew(j,idx)=0; % remove infectiousness for old asymptomatic infection period                     
                       if mig_out && s>rng(j,2)-1 % if j is an internal migrator and recovered during 2nd observation
                           RAobs2actvAnew=RAobs2actvA(RAobs2actvA~=j1); % remove from set of 2nd observations in which individuals are infectious
                           hAnew(j1,idx3)=0; % remove infectiousness from 2nd observation
                           idx=[idx,idx3]; % columns to change in infectiousness matrix
                           tmp=rateHHA(:,ib(j))*hA(j,idx)+rateHHA(:,ib(j1))*hA(j1,idx); % infection pressure from old asymptomatic infection period to subtract
                       else
                           tmp=rateHHA(:,ib(j))*hA(j,idx); % infection pressure to subtract
                       end
                       % Calculate new infection pressure
                       lambdaHH_new(:,idx)=lambdaHH(:,idx)-tmp; % new infection pressure on each HH
                       lambdaHHA_new(:,idx)=lambdaHHA(:,idx)-tmp; % new infection pressure on each HH from asymptomatics
                       lambda_new(:,idx)=lambdaHH_new(ib,idx); % new infection pressure on each individual
                       q=log(prob0(j,2)/sum(prob0(j,2:3)))-log(prob0(j,3)/sum(prob0(j,2:3))); % proposal ratio for recovery time move from after start of study to before
                   end
               else % recovered after start of study
                   if mig_out % if j is internal migrator and this is their 1st observation
                       probAIPnew=[geopdf(0:rng(j1,2)-2,p2),1-geocdf(rng(j1,2)-2,p2)]; % calculate probabilities of different possible infection durations (incl. ones longer than period for which j was present)
                       fwdRA=[1:rng(j1,2)-1,tmax+1]; % possible new recovery times                     
                   else
                       probAIPnew=[geopdf(0:rng(j,2)-2,p2),1-geocdf(rng(j,2)-2,p2)];                       
                       fwdRA=[1:rng(j,2)-1,tmax+1];
                   end
                   sp=fwdRA(sum(rand>=cumsum(probAIPnew))+1); % draw new recovery time
                   if s==0 % if j is not already in asymptomatic infectiousness matrix as they recovered before the start of the study
                      prevAnew=prevA(prevA~=j); % remove them from index of initially previously asymptomatically infected
                      actvAnew=[actvA;j]; % add them to index of initially active asymptomatics
                      A1new=[A1;j]; % add them to index of initially active asymptomatics and individuals asymptomatically infected during the study
                      idx=1:min(min(sp,rng(j,2)),tmax); % index of columns to add infectiousness to in infectiousness matrix
                      hAnew(j,idx)=pold(6); % add infectiousness
                      if mig_out && sp>rng(j,2)-1 % if j is an internal migrator and proposed recovery time is during 2nd observation
                          RAobs2actvAnew=[RAobs2actvA;j1]; % add to set of 2nd observations in which individuals are infectious
                          idx4=rng(j,2)+1:min(min(sp,rng(j1,2)),tmax); % columns to change in infectiousness matrix
                          hAnew(j1,idx4)=pold(6); % add infectiousness to 2nd observation
                          idx=[idx,idx4]; % columns to change in infection pressure matrix
                          tmp=rateHHA(:,ib(j))*hAnew(j,idx)+rateHHA(:,ib(j1))*hAnew(j1,idx); % infection pressure from new asymptomatic infection period to add
                      else
                          tmp=rateHHA(:,ib(j))*hAnew(j,idx);
                      end
                      % Calculate new infection pressure
                      lambdaHH_new(:,idx)=lambdaHH(:,idx)+tmp; % new infection pressure on each HH
                      lambdaHHA_new(:,idx)=lambdaHHA(:,idx)+tmp; % new infection pressure on each HH from asymptomatic 
                      lambda_new(:,idx)=lambdaHH_new(ib,idx); % new infection pressure on each individual
                      q=log(prob0(j,3)/sum(prob0(j,2:3)))-log(prob0(j,2)/sum(prob0(j,2:3))); % proposal ratio for recovery time move from before start of study to after
                   else % j is already in asymptomatic infectiousness matrix as they recovered after start of study
                      idx2=1:min(min(sp,rng(j,2)),tmax); % columns to add infectiousness to in infectiousness matrix
                      hAnew(j,idx1)=0; % remove infectiousness from old infection period
                      hAnew(j,idx2)=pold(6); % add infectiousness from new infection period
                      idx=union(idx1,idx2); % columns to change in infection pressure matrix
                      if ~mig_out || (mig_out && s<=rng(j,2)-1 && sp<=rng(j,2)-1) % if j is not an internal migrator or they are an internal migrator but both old and new recovery times are during 1st observation
                          tmp=rateHHA(:,ib(j))*(hAnew(j,idx)-hA(j,idx)); % infection pressure to add
                      else % j internally migrated and old or new recovery time is after migration time
                          if s<=rng(j,2)-1 && sp>rng(j,2)-1 % old recovery time is after migration, new one is after
                              RAobs2actvAnew=[RAobs2actvA;j1]; % add to index of 2nd observations in which individuals are infectious
                              idx4=rng(j,2)+1:min(min(sp,rng(j1,2)),tmax); % columns to add new infectiousness to in infectiousness matrix
                              hAnew(j1,idx4)=pold(6); % add infectiousness to 2nd observation
                              idx=[idx,idx4]; % columns to change in infection pressure matrix
                              tmp=rateHHA(:,ib(j))*(hAnew(j,idx)-hA(j,idx))+rateHHA(:,ib(j1))*hAnew(j1,idx); % infection pressure to add
                          elseif s>rng(j,2)-1 && sp<=rng(j,2)-1 % old recovery time is after migration, new one is before
                              RAobs2actvAnew=RAobs2actvA(RAobs2actvA~=j1); % add to index of 2nd observations in which individuals are infectious
                              hAnew(j1,idx3)=0; % remove infectiousness if currently still infected after migration
                              idx=[idx,idx3]; % columns to change in infection pressure matrix
                              tmp=rateHHA(:,ib(j))*(hAnew(j,idx)-hA(j,idx))-rateHHA(:,ib(j1))*hA(j1,idx); % infection pressure to add
                          else % both old and new recovery time are after migration
                              hAnew(j1,idx3)=0; % remove infectiousness from 2nd observation if currently still infected after migration
                              idx4=rng(j,2)+1:min(min(sp,rng(j1,2)),tmax); % columns to add infectiousness to in infectiousness matrix
                              hAnew(j1,idx4)=pold(6); % add infectiousness to 2nd observation
                              idx=[idx,union(idx3,idx4)]; % columns to change in infection pressure matrix
                              tmp=rateHHA(:,ib(j))*(hAnew(j,idx)-hA(j,idx))+rateHHA(:,ib(j1))*(hAnew(j1,idx)-hA(j1,idx)); % infection pressure to add
                          end
                      end
                      % Calculate new infection pressure
                      lambdaHH_new(:,idx)=lambdaHH(:,idx)+tmp; % new infection pressure on HH
                      lambdaHHA_new(:,idx)=lambdaHHA(:,idx)+tmp; % new infection pressure on HH from asymptomatic individuals
                      lambda_new(:,idx)=lambdaHH_new(ib,idx); % new infection pressure on individuals
                      q=0; % backward and forward proposal probabilities cancel
                  end              
               end
           elseif tp==tmax+1 % proposed not asymptomatically infected before or during study
               sp=tmax+1; % new recovery time is tmax+1
               bck=[0:rng(j,2)-1,tmax+1]; % possible infection times for reverse move 
               if s==0 % if j is not already in asymptomatic infectiousness matrix as they recovered by start of study
                   prevAnew=prevA(prevA~=j); % remove them from set of initially previously asymptomatically infected
                   idx=[]; % no change required
                   q=log(prob0(j,3)/sum(prob0(j,2:3))); % proposal ratio for recovery time move from before study to no infection (tmax+1)
               else % j is already in asymptomatic infectiousness matrix as they recovered after start of study
                   actvAnew=actvA(actvA~=j); % remove from initially active asymptomatics
                   A1new=A1(A1~=j); % remove them from index of initially active asymptomatics and individuals asymptomatically infected during the study
                   idx=idx1; % columns to remove old infectiousness from in infectiousness matrix
                   hAnew(j,idx)=0; % remove infectiousness
                   if mig_out && s>rng(j,2)-1 % if j is an internal migrator and recovered during 2nd observation
                       RAobs2actvAnew=RAobs2actvA(RAobs2actvA~=j1); % remove from set of 2nd observations in which individuals are infectious
                       hAnew(j1,idx3)=0; % remove infectiousness if currently still infected after migration
                       idx=[idx,idx3]; % columns to change in infection pressure matrix
                       tmp=rateHHA(:,ib(j))*hA(j,idx)+rateHHA(:,ib(j1))*hA(j1,idx); % infection pressure to subtract
                   else
                       tmp=rateHHA(:,ib(j))*hA(j,idx);
                   end
                   % Calculate new infection pressure
                   lambdaHH_new(:,idx)=lambdaHH(:,idx)-tmp; % new infection pressure on each HH
                   lambdaHHA_new(:,idx)=lambdaHHA(:,idx)-tmp; % new infectipn pressure on each HH from asymptomatic individuals
                   lambda_new(:,idx)=lambdaHH_new(ib,idx); % new infection pressure on each individual
                   q=log(prob0(j,2)/sum(prob0(j,2:3))); % proposal ratio for recovery time move from during study to no infection
               end
               S0new=[S0;j]; % add j to initially susceptible set
               if mig_out % if j is an internal migrator and it is their 1st observation 
                   Snew(j1,rng(j,2):rng(j1,2)-1)=1; % add susceptibility to 2nd observation now that individual is not asymptomatically infected before start of study
               end
           else % proposed asymptomatically infected during study (tp in [1,tmax])
               if mig_out % if j is an internal migrator and this is their 1st observation
                   probAIPnew=[geopdf(0:rng(j1,2)-tp-2,p2),1-geocdf(rng(j1,2)-tp-2,p2)]; % probabilities of different infection period durations
                   fwdRA=[tp+1:rng(j1,2)-1,tmax+1]; % possible new recovery times
               else
                   probAIPnew=[geopdf(0:rng(j,2)-tp-2,p2),1-geocdf(rng(j,2)-tp-2,p2)];
                   fwdRA=[tp+1:rng(j,2)-1,tmax+1];
               end
               sp=fwdRA(sum(rand>=cumsum(probAIPnew))+1); % draw new recovery time
               idx2=tp+1:min(min(sp,rng(j,2)),tmax); % columns to add infectiousness to in infectiousness matrix
               if mig_out && sp>rng(j,2)-1 % if j is an internal migrator and new recovery time is during 2nd observation
                   idx4=rng(j,2)+1:min(min(sp,rng(j1,2)),tmax); % columns to add infectiousness for 2nd observation to in infectiousness matrix
               end
               bck=[0,max(1,tp-M):min(rng(j,2)-1,tp+M),tmax+1]; % possible infection times for reverse move
               tAmnew(j,tp)=1; % add new infection time
               if s==0 % if j is not already in asymptomatic infectiousness matrix as they recovered by start of study
                   prevAnew=prevA(prevA~=j); % remove them from previously asymptomatically infected set
                   A1new=[A1;j]; % add them to index of initially active asymptomatics and individuals asymptomatically infected during the study 
                   idx=idx2;
                   hAnew(j,idx)=pold(6); % add infectiousness for new asymptomatic infection period
                   if mig_out && sp>rng(j,2)-1 % if j is an internal migrator and new recovery time is during 2nd observation
                       RAobs2new=[RAobs2;j1]; % add to index of 2nd observations in which individuals are infectious
                       hAnew(j1,idx4)=pold(6); % add infectiousness to 2nd observation 
                       idx=[idx,idx4]; % columns to add infection pressure to in infection pressure matrix
                       tmp=rateHHA(:,ib(j))*hAnew(j,idx)+rateHHA(:,ib(j1))*hAnew(j1,idx); % infection pressure to add
                   else
                       tmp=rateHHA(:,ib(j))*hAnew(j,idx);
                   end
                   lambdaHH_new(:,idx)=lambdaHH(:,idx)+tmp; % new infection pressure on each HH
                   lambdaHHA_new(:,idx)=lambdaHHA(:,idx)+tmp; % new infection pressure on each HH from asymptomatic individuals
                   lambda_new(:,idx)=lambdaHH_new(ib,idx); % new infection pressure on individuals
                   q=log(prob0(j,3)/sum(prob0(j,2:3))); % proposal ratio for recovery time move from before study
               else % individual is already in asymptomatic infectiousness matrix as they recovered after start of study
                   actvAnew=actvA(actvA~=j); % remove from initially actively asymptomatic individuals
                   hAnew(j,idx1)=0; % remove old infectiousness
                   hAnew(j,idx2)=pold(6); % add new infectiousness
                   idx=union(idx1,idx2); % columns to change in infection pressure matrix
                   % Calculate change in infection pressure accounting for
                   % any overlap in infectiousness between observations for
                   % internal migrators
                   if ~mig_out || (mig_out && s<=rng(j,2)-1 && sp<=rng(j,2)-1)
                       tmp=rateHHA(:,ib(j))*(hAnew(j,idx)-hA(j,idx)); 
                   else % individual internally migrated and old or new recovery time is after migration time
                       if s<=rng(j,2)-1 && sp>rng(j,2)-1
                           RAobs2new=[RAobs2;j1];
                           idx4=rng(j,2)+1:min(min(sp,rng(j1,2)),tmax);
                           hAnew(j1,idx4)=pold(6);
                           idx=[idx,idx4];
                           tmp=rateHHA(:,ib(j))*(hAnew(j,idx)-hA(j,idx))+rateHHA(:,ib(j1))*hAnew(j1,idx);
                       elseif s>rng(j,2)-1 && sp<=rng(j,2)-1
                           RAobs2actvAnew=RAobs2actvA(RAobs2actvA~=j1);
                           hAnew(j1,idx3)=0;
                           idx=[idx,idx3];
                           tmp=rateHHA(:,ib(j))*(hAnew(j,idx)-hA(j,idx))-rateHHA(:,ib(j1))*hA(j1,idx);
                       else
                           RAobs2actvAnew=RAobs2actvA(RAobs2actvA~=j1);
                           RAobs2new=[RAobs2;j1];
                           hAnew(j1,idx3)=0; % remove infectiousness if currently still infected after migration
                           idx4=rng(j,2)+1:min(min(sp,rng(j1,2)),tmax);
                           hAnew(j1,idx4)=pold(6);
                           idx=[idx,union(idx3,idx4)];
                           tmp=rateHHA(:,ib(j))*(hAnew(j,idx)-hA(j,idx))+rateHHA(:,ib(j1))*(hAnew(j1,idx)-hA(j1,idx));
                       end
                   end
                   lambdaHH_new(:,idx)=lambdaHH(:,idx)+tmp;
                   lambdaHHA_new(:,idx)=lambdaHHA(:,idx)+tmp;
                   lambda_new(:,idx)=lambdaHH_new(ib,idx);
                   q=log(prob0(j,2)/sum(prob0(j,2:3))); % proposal ratio for recovery time move from after start of study
               end
               S0new=[S0;j];
           end
           % Calculate full proposal ratio for asymptomatic infection period move
           q=q+log(probA(j,1)/sum(probA(j,bck+1)))-log(probA(j,tp+1)/sum(probA(j,fwd+1)));
           
           % Add susceptibility up to new infection time
           Snew(j,1:min(tp-1,rng(j,2)-1))=1;
           
           % Calculate new log-likelihood terms
           LL1new=L1(Snew,lambda_new);
           LL2new=L2(lambda_new,pold(7),tEm);
           LL4new=L4(lambda_new,pold(7),tAmnew);
           LL5new=L5(age,S0new,actvAnew,prevAnew,pold(5),pold(7),p2);
           LLnew=LL1new+LL2new+LL4new+LL5new;
           LLold=LL1old+LL2old+LL4old+LL5old;
           
           % Calculate M-H acceptance probability
           log_ap=LLnew-LLold+q;
           
           if log_ap > log(rand) % accept new values if acceptance probability > rand
               tA(j)=tp; tRA(j)=sp; S(j,:)=Snew(j,:); tAm(j,:)=tAmnew(j,:);
               if mig_out && tp==tmax+1
                   S(j1,:)=Snew(j1,:);
                   tA(j1)=tmax+1; % change asymptomatic infection time for 2nd observation to tmax+1 for susceptibility
                   tRA(j1)=tmax+1; % change asymptomatic recovery time for 2nd observation to tmax+1 for susceptibility
               end
               hA=hAnew; lambda(:,idx)=lambda_new(:,idx); lambdaHH(:,idx)=lambdaHH_new(:,idx);
               lambdaHHA(:,idx)=lambdaHHA_new(:,idx);
               LL1old=LL1new; LL2old=LL2new;
               LL4old=LL4new; LL5old=LL5new;
               A1=A1new; prevA=prevAnew; actvA=actvAnew; S0=S0new;
               RAobs2actvA=RAobs2actvAnew; RAobs2=RAobs2new;
               acc_Arem=acc_Arem+1;
               acc_rem(j)=acc_rem(j)+1;
           else % else reject and keep old values
               Snew(j,:)=S(j,:); tAmnew(j,:)=tAm(j,:);
               if mig_out && tp==tmax+1
                   Snew(j1,:)=S(j1,:);
               end
               hAnew=hA; lambda_new(:,idx)=lambda(:,idx); lambdaHH_new(:,idx)=lambdaHH(:,idx);
               lambdaHHA_new(:,idx)=lambdaHHA(:,idx);
               A1new=A1; prevAnew=prevA; actvAnew=actvA; S0new=S0;
               RAobs2actvAnew=RAobs2actvA; RAobs2new=RAobs2;
               rej_Arem=rej_Arem+1;
               rej_rem(j)=rej_rem(j)+1;
           end
       elseif t==tmax+1 % currently not asymptomatically infected before or during study (N.B. could be 1st or 2nd observation but only matters if it's 1st observation)
           fwd=[rng(j,1):rng(j,2)-1,tmax+1];
           tp=fwd(sum(rand>=cumsum(probA(j,fwd+1))/sum(probA(j,fwd+1)))+1);
           if tp~=tmax+1 % calculate likelihood if proposed infection time is different
              if tp==0 % proposed asymptomatically infected before start of study (tp=0) (N.B. must be 1st observation)
                  bck=[0:rng(j,2)-1,tmax+1];                  
                  if rand<prob0(j,3)/sum(prob0(j,2:3)) % recovered before start of study 
                      sp=0;
                      prevAnew=[prevA;j]; % add j to previously asymptomatically infected set
                      idx=[];
                      q=-log(prob0(j,3)/sum(prob0(j,2:3)));
                  else % recovered after start of study
                      if mig_out
                          probAIPnew=[geopdf(0:rng(j1,2)-2,p2),1-geocdf(rng(j1,2)-2,p2)];
                          fwdRA=[1:rng(j1,2)-1,tmax+1];
                      else
                          probAIPnew=[geopdf(0:rng(j,2)-2,p2),1-geocdf(rng(j,2)-2,p2)];
                          fwdRA=[1:rng(j,2)-1,tmax+1];
                      end
                      sp=fwdRA(sum(rand>=cumsum(probAIPnew))+1);
                      idx=1:min(min(sp,rng(j,2)),tmax);
                      A1new=[A1;j];
                      actvAnew=[actvA;j];
                      hAnew(j,idx)=pold(6);
                      if mig_out && sp>rng(j,2)-1
                          RAobs2actvAnew=[RAobs2actvA;j1];
                          idx4=rng(j,2)+1:min(min(sp,rng(j1,2)),tmax);
                          hAnew(j1,idx4)=pold(6);
                          idx=[idx,idx4];
                          tmp=rateHHA(:,ib(j))*hAnew(j,idx)+rateHHA(:,ib(j1))*hAnew(j1,idx);
                      else
                          tmp=rateHHA(:,ib(j))*hAnew(j,idx);
                      end                      
                      lambdaHH_new(:,idx)=lambdaHH(:,idx)+tmp;
                      lambdaHHA_new(:,idx)=lambdaHHA(:,idx)+tmp;
                      lambda_new(:,idx)=lambdaHH_new(ib,idx);
                      q=-log(prob0(j,2)/sum(prob0(j,2:3)));
                  end
                  S0new=S0(S0~=j); % remove j from initially susceptible set
              else % proposed asymptomatically infected during study (tp in [1,tmax])
                  if mig_out
                      probAIPnew=[geopdf(0:rng(j1,2)-tp-2,p2),1-geocdf(rng(j1,2)-tp-2,p2)];
                      fwdRA=[tp+1:rng(j1,2)-1,tmax+1];
                  else
                      probAIPnew=[geopdf(0:rng(j,2)-tp-2,p2),1-geocdf(rng(j,2)-tp-2,p2)];
                      fwdRA=[tp+1:rng(j,2)-1,tmax+1];
                  end                  
                  sp=fwdRA(sum(rand>=cumsum(probAIPnew))+1);
                  idx=tp+1:min(min(sp,rng(j,2)),tmax);
                  if rng(j,1)==0 % asymptomatic infection before study is possible as individual was alive/present
                      bck=[0,max(1,tp-M):min(rng(j,2)-1,tp+M),tmax+1];
                  else % individual was born or imigrated after start of study
                      bck=[max(rng(j,1),tp-M):min(rng(j,2)-1,tp+M),tmax+1];
                  end
                  tAmnew(j,tp)=1; % add new infection time
                  A1new=[A1;j];
                  hAnew(j,idx)=pold(6);
                  if mig_out && sp>rng(j,2)-1
                      RAobs2new=[RAobs2;j1];
                      idx4=rng(j,2)+1:min(min(sp,rng(j1,2)),tmax);
                      hAnew(j1,idx4)=pold(6);
                      idx=[idx,idx4];
                      tmp=rateHHA(:,ib(j))*hAnew(j,idx)+rateHHA(:,ib(j1))*hAnew(j1,idx);
                  else
                      tmp=rateHHA(:,ib(j))*hAnew(j,idx);
                  end
                  lambdaHH_new(:,idx)=lambdaHH(:,idx)+tmp;
                  lambdaHHA_new(:,idx)=lambdaHHA(:,idx)+tmp;
                  lambda_new(:,idx)=lambdaHH_new(ib,idx);
                  q=0;
              end
              q=q+log(probA(j,tmax+2)/sum(probA(j,bck+1)))-log(probA(j,tp+1)/sum(probA(j,fwd+1))); % calculate proposal ratio

              Snew(j,max(1,tp):rng(j,2)-1)=0; % remove susceptibility from new infection time onwards
              if mig_out
                  Snew(j1,rng(j,2):rng(j1,2)-1)=0; % remove susceptibility from 2nd observation now that individual is asymptomatically infected during 1st observation
              end

              LL1new=L1(Snew,lambda_new);
              LL2new=L2(lambda_new,pold(7),tEm);
              LL4new=L4(lambda_new,pold(7),tAmnew);
              LL5new=L5(age,S0new,actvAnew,prevAnew,pold(5),pold(7),p2);
              LLnew=LL1new+LL2new+LL4new+LL5new;
              LLold=LL1old+LL2old+LL4old+LL5old;

              log_ap=LLnew-LLold+q; % calculate M-H acceptance probability

              if log_ap > log(rand)
                  tA(j)=tp; tRA(j)=sp; S(j,:)=Snew(j,:); tAm(j,:)=tAmnew(j,:);
                  if mig_out
                      S(j1,:)=Snew(j1,:);
                      tA(j1)=tmax+2; % change asymptomatic infection time to dummy time for asymptomatic infection not possible now that individual is asymptomatically infected during 1st observation
                      tRA(j1)=tmax+2; % change asymptomatic recovery time to dummy time for asymptomatic infection not possible now that individual is asymptomatically infected during 1st observation
                  end
                  hA=hAnew; lambda(:,idx)=lambda_new(:,idx); lambdaHH(:,idx)=lambdaHH_new(:,idx);
                  lambdaHHA(:,idx)=lambdaHHA_new(:,idx);
                  LL1old=LL1new; LL2old=LL2new;
                  LL4old=LL4new; LL5old=LL5new;
                  A1=A1new; prevA=prevAnew; actvA=actvAnew; S0=S0new;
                  RAobs2actvA=RAobs2actvAnew; RAobs2=RAobs2new;
                  acc_Aadd=acc_Aadd+1;
                  acc_add(j)=acc_add(j)+1;
              else
                  Snew(j,:)=S(j,:); tAmnew(j,:)=tAm(j,:);
                  if mig_out
                      Snew(j1,:)=S(j1,:);
                  end
                  hAnew=hA; lambda_new(:,idx)=lambda(:,idx); lambdaHH_new(:,idx)=lambdaHH(:,idx);
                  lambdaHHA_new(:,idx)=lambdaHHA(:,idx);
                  A1new=A1; prevAnew=prevA; actvAnew=actvA; S0new=S0;
                  RAobs2actvAnew=RAobs2actvA; RAobs2new=RAobs2;
                  rej_Aadd=rej_Aadd+1;
                  rej_add(j)=rej_add(j)+1;
              end
           else % otherwise accept immediately as likelihood doesn't change (if tp=tmax+1)
               acc_Aadd=acc_Aadd+1;
           end
       else % currently asymptomatically infected between time 1 and tmax (N.B. could be 1st or 2nd observation and need to account for both)
           if rng(j,1)==0
               fwd=[0,max(1,t-M):min(rng(j,2)-1,t+M),tmax+1];
           else
               fwd=[max(rng(j,1),t-M):min(rng(j,2)-1,t+M),tmax+1];
           end
           tp=fwd(sum(rand>=cumsum(probA(j,fwd+1))/sum(probA(j,fwd+1)))+1);
           tAmnew(j,t)=0; % remove old infection time
           hAnew(j,idx1)=0; % remove old infectiousness
           if tp==0 % proposed asymptomatically infected before study (N.B. must be 1st observation)
               bck=[0:rng(j,2)-1,tmax+1];
               S0new=S0(S0~=j); % remove j from initially susceptible set
               if rand<prob0(j,3)/sum(prob0(j,2:3)) % recovered before start of study
                   sp=0;
                   A1new=A1(A1~=j);
                   prevAnew=[prevA;j]; % add j to previously asymptomatically infected set                   
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
                   q=-log(prob0(j,3)/sum(prob0(j,2:3)));
               else % recovered after start of study
                   if mig_out
                       probAIPnew=[geopdf(0:rng(j1,2)-2,p2),1-geocdf(rng(j1,2)-2,p2)];
                       fwdRA=[1:rng(j1,2)-1,tmax+1];
                   else
                       probAIPnew=[geopdf(0:rng(j,2)-2,p2),1-geocdf(rng(j,2)-2,p2)];
                       fwdRA=[1:rng(j,2)-1,tmax+1];
                   end
                   sp=fwdRA(sum(rand>=cumsum(probAIPnew))+1);
                   idx2=1:min(min(sp,rng(j,2)),tmax);
                   actvAnew=[actvA;j];
                   hAnew(j,idx2)=pold(6); % add new infectiousness
                   idx=union(idx1,idx2); % create union of indices for removed infectiousness and added infectiousness
                   if ~mig_out || (mig_out && s<=rng(j,2)-1 && sp<=rng(j,2)-1)
                       tmp=rateHHA(:,ib(j))*(hAnew(j,idx)-hA(j,idx));
                   else
                       if s<=rng(j,2)-1 && sp>rng(j,2)-1
                           RAobs2actvAnew=[RAobs2actvA;j1];
                           idx4=rng(j,2)+1:min(min(sp,rng(j1,2)),tmax);
                           hAnew(j1,idx4)=pold(6);
                           idx=[idx,idx4];
                           tmp=rateHHA(:,ib(j))*(hAnew(j,idx)-hA(j,idx))+rateHHA(:,ib(j1))*hAnew(j1,idx);
                       elseif s>rng(j,2)-1 && sp<=rng(j,2)-1
                           RAobs2new=RAobs2(RAobs2~=j1);
                           hAnew(j1,idx3)=0;
                           idx=[idx,idx3];
                           tmp=rateHHA(:,ib(j))*(hAnew(j,idx)-hA(j,idx))-rateHHA(:,ib(j1))*hA(j1,idx);
                       else
                           RAobs2new=RAobs2(RAobs2~=j1);
                           RAobs2actvAnew=[RAobs2actvA;j1];
                           hAnew(j1,idx3)=0;
                           idx4=rng(j,2)+1:min(min(sp,rng(j1,2)),tmax);
                           hAnew(j1,idx4)=pold(6);
                           idx=[idx,union(idx3,idx4)];
                           tmp=rateHHA(:,ib(j))*(hAnew(j,idx)-hA(j,idx))+rateHHA(:,ib(j1))*(hAnew(j1,idx)-hA(j1,idx));
                       end
                   end
                   lambdaHH_new(:,idx)=lambdaHH(:,idx)+tmp;
                   lambdaHHA_new(:,idx)=lambdaHHA(:,idx)+tmp;
                   lambda_new(:,idx)=lambdaHH_new(ib,idx);
                   q=-log(prob0(j,2)/sum(prob0(j,2:3)));
               end
           elseif tp==tmax+1 % proposed not asymptomatically infected before or during study
               sp=tmax+1;
               bck=[rng(j,1):rng(j,2)-1,tmax+1];
               A1new=A1(A1~=j);
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
               q=0;
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
               sp=fwdRA(sum(rand>=cumsum(probAIPnew))+1);
               idx2=tp+1:min(min(sp,rng(j,2)),tmax);
               if rng(j,1)==0 % asymptomatic infection before study is possible as individual was alive/present
                   bck=[0,max(1,tp-M):min(rng(j,2)-1,tp+M),tmax+1];
               else % proposed asymptomatically infected during study and individual was born or imigrated after start of study
                   bck=[max(rng(j,1),tp-M):min(rng(j,2)-1,tp+M),tmax+1];
               end
               tAmnew(j,tp)=1; % add new infection time        
               hAnew(j,idx2)=pold(6);
               idx=union(idx1,idx2);
               if ~mig_out || (mig_out && s<=rng(j,2)-1 && sp<=rng(j,2)-1)
                   tmp=rateHHA(:,ib(j))*(hAnew(j,idx)-hA(j,idx));
               else
                   if s<=rng(j,2)-1 && sp>rng(j,2)-1
                       RAobs2new=[RAobs2;j1];
                       idx4=rng(j,2)+1:min(min(sp,rng(j1,2)),tmax);
                       hAnew(j1,idx4)=pold(6);
                       idx=[idx,idx4];
                       tmp=rateHHA(:,ib(j))*(hAnew(j,idx)-hA(j,idx))+rateHHA(:,ib(j1))*hAnew(j1,idx);
                   elseif s>rng(j,2)-1 && sp<=rng(j,2)-1
                       RAobs2new=RAobs2(RAobs2~=j1);
                       hAnew(j1,idx3)=0;
                       idx=[idx,idx3];
                       tmp=rateHHA(:,ib(j))*(hAnew(j,idx)-hA(j,idx))-rateHHA(:,ib(j1))*hA(j1,idx);
                   else % if s>rng(j,2)-1 && sp>rng(j,2)-1
                       hAnew(j1,idx3)=0;
                       idx4=rng(j,2)+1:min(min(sp,rng(j1,2)),tmax);
                       hAnew(j1,idx4)=pold(6);
                       idx=[idx,union(idx3,idx4)];
                       tmp=rateHHA(:,ib(j))*(hAnew(j,idx)-hA(j,idx))+rateHHA(:,ib(j1))*(hAnew(j1,idx)-hA(j1,idx));
                   end
               end
               lambdaHH_new(:,idx)=lambdaHH(:,idx)+tmp;
               lambdaHHA_new(:,idx)=lambdaHHA(:,idx)+tmp;
               lambda_new(:,idx)=lambdaHH_new(ib,idx);
               q=0;
           end
           q=q+log(probA(j,t+1)/sum(probA(j,bck+1)))-log(probA(j,tp+1)/sum(probA(j,fwd+1))); % calculate proposal ratio
           
           Snew(j,max(1,tp):t-1)=0; % remove susceptibility if new infection time is earlier
           Snew(j,t:min(tp-1,rng(j,2)-1))=1; % add susceptibility if new infection time is later
           
           LL1new=L1(Snew,lambda_new);
           LL2new=L2(lambda_new,pold(7),tEm);
           LL4new=L4(lambda_new,pold(7),tAmnew);
           LL5new=L5(age,S0new,actvAnew,prevAnew,pold(5),pold(7),p2);
           LLnew=LL1new+LL2new+LL4new+LL5new;
           LLold=LL1old+LL2old+LL4old+LL5old;
           
           log_ap=LLnew-LLold+q; % calculate M-H acceptance probability
           
           if log_ap > log(rand)
               tA(j)=tp; tRA(j)=sp; S(j,:)=Snew(j,:); tAm(j,:)=tAmnew(j,:);
               if mig_out && tp==tmax+1
                   S(j1,:)=Snew(j1,:);
                   tA(j1)=tmax+1;
                   tRA(j1)=tmax+1;
               end
               hA=hAnew; lambda(:,idx)=lambda_new(:,idx); lambdaHH(:,idx)=lambdaHH_new(:,idx);
               lambdaHHA(:,idx)=lambdaHHA_new(:,idx);
               LL1old=LL1new; LL2old=LL2new;
               LL4old=LL4new; LL5old=LL5new;
               A1=A1new; prevA=prevAnew; actvA=actvAnew; S0=S0new;
               RAobs2actvA=RAobs2actvAnew; RAobs2=RAobs2new;
               acc_Amov=acc_Amov+1;
               acc_mov(j)=acc_mov(j)+1;
           else
               Snew(j,:)=S(j,:); tAmnew(j,:)=tAm(j,:);
               if mig_out && tp==tmax+1
                   Snew(j1,:)=S(j1,:);
               end
               hAnew=hA; lambda_new(:,idx)=lambda(:,idx); lambdaHH_new(:,idx)=lambdaHH(:,idx);
               lambdaHHA_new(:,idx)=lambdaHHA(:,idx);
               A1new=A1; prevAnew=prevA; actvAnew=actvA; S0new=S0;
               RAobs2actvAnew=RAobs2actvA; RAobs2new=RAobs2;
               rej_Amov=rej_Amov+1;
               rej_mov(j)=rej_mov(j)+1;
           end
       end
       end
    end

    %% UPDATE ASYMPTOMATIC INFECTION AND RECOVERY TIMES OF PKDL CASES WITHOUT PRIOR KA
    for i=1:nPA
        j=PA(i); % get index of PKDL case without prior KA
        t=tA(j); % get current asymptomatic infection time of j
        s=tRA(j); % get current asymptomatic recovery time of j
%         probDIP=[nbinpdf(0:tP(j)-1,r3,p3),1-nbincdf(tP(j)-1,r3,p3)];
%         sp=tP(j)-randsample(tP(j),1,true,probDIP);
        sp=tP(j)-nbinrnd(r3,p3);
        tp=sp-(geornd(p2)+1);
%         if sp<0
%             sp=0;
%         end
        
        if tp>=max(1,tB(j)+1) % propose new infection time if asymptomatically infected after start of study and birth month
%             if tRA(j)>0 || sp>0
%                 if t>0 % currently asymptomatically infected during study
                    tAmnew(j,t)=0; % remove old infection time if it is after start of study
%                     if tp<=0 && sp>0 % proposed initially actively asymptomatically infected
%                         actvAPAnew=[actvAPA;j];
%                         S0PAnew=S0PA(S0PA~=j);
%                     elseif sp<=0 % proposed recovered before study
%                         prevAPAnew=[prevAPA;j];
%                         S0PAnew=S0PA(S0PA~=j);
%                     end
%                 else % currently asymptomatically infected before study (t<=0)
%                     if tRA(j)>0 % currently initially actively asymptomatically infected
%                         if tp>0 % proposed asymptomatically infected during study
%                             actvAPAnew=actvAPA(actvAPA~=j);
%                             S0PAnew=[S0PA;j];    
%                         elseif sp<=0 % proposed asymptomatically infected before study
%                             prevAPAnew=[prevAPA;j];
%                             actvAPAnew=actvAPA(actvAPA~=j);
%                         end
%                     else % tRA(j)<=0
%                         if tp<=0 && sp>0
%                             prevAPAnew=prevAPA(prevAPA~=j);
%                             actvAPAnew=[actvAPA;j];
%                         elseif tp>0
%                             prevAPAnew=prevAPA(prevAPA~=j);
%                             S0PAnew=[S0PA;j];
%                         end
%                     end
%                 end
%                 if tp>0
%                     tAmnew(j,tp)=1; % add new infection time if it is after start of study
%                 end
                % Update infection time and susceptibility matrices
                tAmnew(j,tp)=1; % add new infection time to infection time matrix
                Snew(j,max(1,tp):t-1)=0; % remove old susceptible times if new infection time is earlier (from month indvdl is infected up to but not incl. month before old infection time)
                Snew(j,max(1,t):tp-1)=1; % add new susceptible times if new infection time is later (up to but not incl. month individual is infected)
                
                % Update infectiousness matrix
                hPAnew(i,max(0,t)+1:s)=0; % remove old infectiousness
                hPAnew(i,max(0,tp)+1:sp)=h40; % add new infectiousness
                erlrA=max(0,min(t,tp)); % index of earlier infection time between old and new infection time
                ltrA=max(0,max(t,tp)); % index of later infection time between old and new infection time
                erlrRA=max(0,min(tRA(j),sp)); % index of earlier recovery time between old and new recovery time
                ltrRA=max(tRA(j),sp); % index of later recovery time between old and new recovery time
                
                % Update infection pressure
                idx=[erlrA+1:ltrA,erlrRA+1:ltrRA]; % index of columns to change
                lambdaHHPA_new(:,idx)=lambdaHHPA(:,idx)-rateHHPA(:,i)*hPA(i,idx)+rateHHPA(:,i)*hPAnew(i,idx); % new infection pressure from PKDL cases w/o prior KA
                lambdaHH_new(:,idx)=lambdaHH(:,idx)-rateHHPA(:,i)*hPA(i,idx)+rateHHPA(:,i)*hPAnew(i,idx); % new infection pressure on each HH
                lambda_new(:,idx)=lambdaHH_new(ib,idx); % new infection pressure on each individual
                
                % Calculate new log-likelihood terms
                LL1new=L1(Snew,lambda_new);
                LL2new=L2(lambda_new,pold(7),tEm);
                LL4new=L4(lambda_new,pold(7),tAmnew);
                LL6new=L5(age,S0PAnew,actvAPAnew,prevAPAnew,pold(5),pold(7),p2);
                LLnew=LL1new+LL2new+LL4new+LL6new;
                LLold=LL1old+LL2old+LL4old+LL6old;
                
                % Calculate Metropolis-Hastings acceptance probability
                log_ap=LLnew-LLold;
            
                if log_ap > log(rand) % accept new values if acceptance probability > rand
                    tA(j)=tp; tRA(j)=sp; S(j,:)=Snew(j,:); tAm(j,:)=tAmnew(j,:);
                    hPA=hPAnew; lambda(:,idx)=lambda_new(:,idx); lambdaHH(:,idx)=lambdaHH_new(:,idx);
                    lambdaHHPA(:,idx)=lambdaHHPA_new(:,idx);
                    LL1old=LL1new; LL2old=LL2new;
                    LL4old=LL4new; LL6old=LL6new;
                    prevAPA=prevAPAnew; S0PA=S0PAnew; actvAPA=actvAPAnew;
                    acc_PA=acc_PA+1;
                else % else reject and keep old values
                    Snew(j,:)=S(j,:); tAmnew(j,:)=tAm(j,:);
                    hPAnew=hPA; lambda_new(:,idx)=lambda(:,idx); lambdaHH_new(:,idx)=lambdaHH(:,idx);
                    lambdaHHPA_new(:,idx)=lambdaHHPA(:,idx);
                    prevAPAnew=prevAPA; S0PAnew=S0PA; actvAPAnew=actvAPA;
                    rej_PA=rej_PA+1;
                end
%             else % if both old and new recovery times are in/before start month then likelihood does not change so always accept
%                 tA(j)=tp; tRA(j)=sp;
%                 acc_PA=acc_PA+1;
%             end
        else % otherwise reject immediately
            rej_PA=rej_PA+1;
        end
    end

    %% UPDATE MISSING KA ONSET TIMES
    for i=1:nNO
        j=NO(i); % get index of KA case without onset time
        tInew(j)=tR(j)-(nbinrnd(r0,p0)+1); % draw new onset time by drawing new onset-to-treatment time and subtracting it from the treatment time
        
        if tInew(j)>tE(j) && tInew(j)>=tIlb(j) && tInew(j)<=tIub(j) % calculate log-likelihood if new onset time is after infection time, within infection time bounds, and before death
            tIj=tI(j);
            tIjnew=tInew(j);
            IPnew(j)=tIjnew-tE(j); % new incubation period
            m=(IPNIA==j); % get index of case in infectiousness matrix
            hnew(m,tIj+1:tIjnew)=h0; % reduce infectiousness up to new onset time if it is later
            hnew(m,tIjnew+1:tIj)=1; % increase infectiousness from new onset time if it is earlier
            erlrI=min(tIj,tIjnew); % index of column for earlier onset time between old and new onset time
            ltrI=max(tIj,tIjnew); % index of column for later onset time between old and new onset time
            
            % Calculate new infection pressure
            idx=erlrI+1:ltrI; % index of columns to change
            lambdaHH_new(:,idx)=rateHH*hnew(:,idx)+lambdaHHPA(:,idx)+lambdaHHA(:,idx)+pold(3); % new infection pressure on each HH
            lambda_new(:,idx)=lambdaHH_new(ib,idx); % infection pressure on each individual

            % Calculate new log-likelihood terms
            LL1new=L1(S,lambda_new);
            LL2new=L2(lambda_new,pold(7),tEm);
            LL3new=L3(IPnew(I),r1,p1new);
            LL4new=L4(lambda_new,pold(7),tAm);
            LLnew=LL1new+LL2new+LL3new+LL4new;
            LLold=LL1old+LL2old+LL3old+LL4old;
            
            % Calculate Metropolis-Hastings acceptance probability
            log_ap=LLnew-LLold;
            
            if log_ap > log(rand) % accept new values if acceptance probability > rand
                IPold(j)=IPnew(j); tI(j)=tInew(j); h(m,:)=hnew(m,:); lambda(:,idx)=lambda_new(:,idx); lambdaHH(:,idx)=lambdaHH_new(:,idx);
                LL1old=LL1new; LL2old=LL2new; LL3old=LL3new; 
                LL4old=LL4new; % keep updated info
                acc_I=acc_I+1;
            else % else reject and keep new values
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
        j=NONR(i); % get index of KA case without onset or treatment times
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
            tIj=tI(j);
            tIjnew=tInew(j);
            tRj=tR(j);
            tRjnew=tRnew(j);   
            
            % calculate log-likelihood if an infection-treatment block move 
            % is possible; allow treatment month to be same as death month 
            % as cases were said to have died during treatment
            if tIj>maxIP && tIjnew>maxIP % only update infection time and susceptibility matrices for cases with onset after initial window and after entry
                tEmnew(j,tEj)=0; % remove old infection time
                tEmnew(j,tEjnew)=1; % add new infection time
                Snew(j,tEjnew:tEj-1)=0; % remove old susceptible times if new infection time is earlier (from month individual is infected up to but not incl. month before old infection time)
                Snew(j,tEj:tEjnew-1)=1; % add new susceptible times if new infection time is later (up to but not incl. month individual is infected)
            elseif tIj>maxIP && tIjnew<=maxIP % if old onset time is after incubation period window and new one is within it
                tEmnew(j,tEj)=0; % remove infection time from infection time matrix
                Snew(j,:)=0; % remove susceptibility for case
            elseif tIj<=maxIP && tIjnew>maxIP % if old onset time was within incubation period window and new one is after it
                tEmnew(j,tEjnew)=1; % add new infection time to infection time matrix
                Snew(j,1:tEjnew-1)=1; % add susceptibility up to new infection time
            end
                     
            m=(IPNIA==j); % get index of case in infectiousness matrix
            hnew(m,max(0,tEj)+1:tIj)=0; % remove old presymptomatic infectiousness
            hnew(m,tIj+1:tRj)=0; % remove old symptomatic infectiousness
            hnew(m,max(0,tEjnew)+1:tIjnew)=h0; % add new presymptomatic infectiousness
            hnew(m,tIjnew+1:tRjnew)=1; % add new symptomatic infectiousness
  
            erlrE=max(0,min(tEj,tEjnew)); % index of column for earlier infection time between old and new infection time
            ltrE=max(0,max(tEj,tEjnew)); % index of column for later infection time between old and new infection time
            erlrI=min(tIj,tIjnew); % index of column for earlier onset time between old and new onset time
            ltrI=max(tIj,tIjnew); % index of column for later onset time between old and new onset time
            erlrR=min(tRj,tRjnew); % index of column for earlier treatment time between old and new treatment time
            ltrR=max(tRj,tRjnew); % index of column for later treatment time between old and new treatment time
            
            % Calculate new infection pressure
            idx=[erlrE+1:ltrE,erlrI+1:ltrI,erlrR+1:ltrR]; % index of columns to change
            lambdaHH_new(:,idx)=rateHH*hnew(:,idx)+lambdaHHPA(:,idx)+lambdaHHA(:,idx)+pold(3); % new infection pressure on each HH
            lambda_new(:,idx)=lambdaHH_new(ib,idx); % new infection pressure on each individual
            
            % Calculate new log-likelihood terms
            LL1new=L1(Snew,lambda_new);
            LL2new=L2(lambda_new,pold(7),tEmnew);
            LL4new=L4(lambda_new,pold(7),tAm);
            LLnew=LL1new+LL2new+LL4new;
            LLold=LL1old+LL2old+LL4old;
            
            % Calculate Metropolis-Hastings acceptance probability
            log_ap=LLnew-LLold;
            
            if log_ap > log(rand) % accept new values if acceptance probability > rand
                tE(j)=tEnew(j); S(j,:)=Snew(j,:); tEm(j,:)=tEmnew(j,:); tI(j)=tInew(j); tR(j)=tRnew(j);
                h(m,:)=hnew(m,:); lambda(:,idx)=lambda_new(:,idx); lambdaHH(:,idx)=lambdaHH_new(:,idx);
                LL1old=LL1new; LL2old=LL2new; 
                LL4old=LL4new; % keep updated info
                acc_ERmove=acc_ERmove+1;
            else % else reject and keep old values
                tEnew(j)=tE(j); Snew(j,:)=S(j,:); tEmnew(j,:)=tEm(j,:); tInew(j)=tI(j); tRnew(j)=tR(j);
                hnew(m,:)=h(m,:); lambda_new(:,idx)=lambda(:,idx); lambdaHH_new(:,idx)=lambdaHH(:,idx); % keep old values, don't change log-likelihood
                rej_ERmove=rej_ERmove+1;
            end
        else % otherwise reject immediately
            tEnew(j)=tE(j); tInew(j)=tI(j); tRnew(j)=tR(j);
            rej_ERmove=rej_ERmove+1;
        end
    end
    
    %% UPDATE MISSING KA TREATMENT TIMES
    for i=1:nONR
        j=ONR(i); % get index of case without treatment time
        tRnew(j)=tI(j)+nbinrnd(r0,p0)+1; % draw new treatment time by drawing new onset-to-treatment time and adding it to onset time
              
        if tRnew(j)<=min(min(tP(j)-1,tD(j)),tmax) % calculate log-likelihood if new recovery time is before PKDL onset/death/end of study
            tRj=tR(j);
            tRjnew=tRnew(j);
            m=(IPNIA==j); % get infex of case in infectiousness matrix
            hnew(m,tRj+1:tRjnew)=1; % add infectiousness up to new treatment time if it is later
            hnew(m,tRjnew+1:tRj)=0; % remove infectiousness from new treatment time if it is earlier
            erlrR=min(tRj,tRjnew); % index of column for earlier treatment time between old and new treatment time
            ltrR=max(tRj,tRjnew); % index of column for later treatment time between old and new treatment time
            
            % Calculate new infection pressure
            idx=erlrR+1:ltrR; % index of columns to change
            lambdaHH_new(:,idx)=rateHH*hnew(:,idx)+lambdaHHPA(:,idx)+lambdaHHA(:,idx)+pold(3); % new infection pressure on each HH
            lambda_new(:,idx)=lambdaHH_new(ib,idx); % new infection pressure on each individual
            
            % Calculate new log-likelihood terms
            LL1new=L1(S,lambda_new);
            LL2new=L2(lambda_new,pold(7),tEm);
            LL4new=L4(lambda_new,pold(7),tAm);
            LLnew=LL1new+LL2new+LL4new;
            LLold=LL1old+LL2old+LL4old;
            
            % Calculate Metropolis-Hastings acceptance value
            log_ap=LLnew-LLold;
            
            if log_ap > log(rand) % accept new values if acceptance probability > rand
                tR(j)=tRnew(j); h(m,:)=hnew(m,:); lambda(:,idx)=lambda_new(:,idx); lambdaHH(:,idx)=lambdaHH_new(:,idx);
                LL1old=LL1new; LL2old=LL2new; 
                LL4old=LL4new; % keep updated info
                acc_R=acc_R+1;
            else % else reject and keep old values
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
        j=ANONR(i); % get index of potentially initially active KA case
        tmp2=0; % onset-to-treatment period shift
        while tmp2==0 % resample until a non-zero value is obtained
            tmp2=round(sqrt(ERvar)*randn);
        end
        tInew(j)=tI(j)+tmp2; % new onset time
        tRnew(j)=tInew(j)+nbinrnd(r0,p0)+1; % draw new treatment time by drawing new onset-to-treatment period and adding it to new onset time
            
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
                m=(IPNIA==j); % get index of case in infectiousness matrix
                hnew(m,max(0,tRj)+1:tRjnew)=1; % add infectiousness up to new treatment time if it is later
                hnew(m,max(0,tRjnew)+1:tRj)=0; % remove infectiousness from new treatment time if it is later
                erlrR=max(0,min(tRj,tRjnew)); % index of column for earlier treatment time between old and new treatment time
                ltrR=max(tRj,tRjnew); % index of column for later treatment time between old and new treatment time
                
                % Calculate new infection pressure
                idx=erlrR+1:ltrR; % index of columns to change
                lambdaHH_new(:,idx)=rateHH*hnew(:,idx)+lambdaHHPA(:,idx)+lambdaHHA(:,idx)+pold(3); % new infection pressure on each HH
                lambda_new(:,idx)=lambdaHH_new(ib,idx); % new infection pressure on each individual
                
                % Calculate new log-likelihood terms
                LL1new=L1(S,lambda_new);
                LL2new=L2(lambda_new,pold(7),tEm);
                LL4new=L4(lambda_new,pold(7),tAm);
                LLnew=LL1new+LL2new+LL4new;
                LLold=LL1old+LL2old+LL4old;
                
                % Calculate Metropolis-Hastings acceptance probability
                log_ap=LLnew-LLold;
                
                if log_ap > log(rand) % accept new values if acceptance probability > rand
                    tI(j)=tInew(j); tR(j)=tRnew(j);
                    h(m,:)=hnew(m,:); lambda(:,idx)=lambda_new(:,idx); lambdaHH(:,idx)=lambdaHH_new(:,idx); 
                    LL1old=LL1new; LL2old=LL2new; 
                    LL4old=LL4new; % keep updated info
                    acc_AIRmove=acc_AIRmove+1;
                else % else reject and keep old values
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
        j=AONR(i); % get index of potentially initially active KA case
        tRnew(j)=tI(j)+nbinrnd(r0,p0)+1; % draw new treatment time by drawing new onset-to-treatment time and adding it to onset time
        
        if tRnew(j)<=min(min(tP(j)-1,tD(j)),tmax) % calculate log-likelihood if new recovery time is before PKDL onset/death/end of study
            tRj=tR(j);
            tRjnew=tRnew(j);           
            if tRj>0 || tRjnew>0 % calculate log-likelihood if old or new recovery time is after start month
                m=(IPNIA==j); % get index of case in infectiousness matrix
                hnew(m,max(0,tRj)+1:tRjnew)=1; % add infectiousness up to new treatment time if it is later
                hnew(m,max(0,tRjnew)+1:tRj)=0; % remove infectiousness from new treatment time if it is earlier
                erlrR=max(0,min(tRj,tRjnew)); % index of column for earlier treatment time between old and new treatment time
                ltrR=max(tRj,tRjnew); % index of column for later treatment time between old and new treatment time
                
                % Calculate new infection pressure
                idx=erlrR+1:ltrR; % index of columns to change
                lambdaHH_new(:,idx)=rateHH*hnew(:,idx)+lambdaHHPA(:,idx)+lambdaHHA(:,idx)+pold(3); % new infection pressure on each HH
                lambda_new(:,idx)=lambdaHH_new(ib,idx); % new infection pressure on each individual
                
                % Calculate new log-likelihood terms
                LL1new=L1(S,lambda_new);
                LL2new=L2(lambda_new,pold(7),tEm);
                LL4new=L4(lambda_new,pold(7),tAm);
                LLnew=LL1new+LL2new+LL4new;
                LLold=LL1old+LL2old+LL4old;                
                
                % Calculate Metropolis-Hastings acceptance value
                log_ap=LLnew-LLold;
                
                if log_ap > log(rand) % accept new values if acceptance probability > rand
                    tR(j)=tRnew(j); h(m,:)=hnew(m,:); lambda(:,idx)=lambda_new(:,idx); lambdaHH(:,idx)=lambdaHH_new(:,idx); 
                    LL1old=LL1new; LL2old=LL2new;
                    LL4old=LL4new; % keep updated info
                    acc_AR=acc_AR+1;
                else % else reject and keep old values
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
       j=RLNO(i); % get index of relapse case without relapse or relapse treatment times
       tRLnew(j)=tR(j)+geornd(p4)+1; % draw new relapse time by drawing new treatment-to-relapse time and adding it to treatment time
       tRLRnew(j)=tRLnew(j)+nbinrnd(r0,p0)+1; % draw new relapse treatment time by drawing time from KA onset-to-treatment time distribution and adding it to relapse time 
       
       if tRLnew(j)<=min(min(min(tEM(j),tP(j)),tD(j))-2,tmax) && tRLRnew(j)<=min(min(min(tEM(j),tP(j)),tD(j))-1,tmax)
           tRLj=tRL(j);
           tRLjnew=tRLnew(j);
           tRLRj=tRLR(j);
           tRLRjnew=tRLRnew(j);
           m=(IPNIA==j); % get index of case in infectiousness matrix
           hnew(m,tRLj+1:tRLRj)=0; % remove old infectiousness
           hnew(m,tRLjnew+1:tRLRjnew)=1; % add new infectiousness
           erlrRL=min(tRLj,tRLjnew); % index of column for earlier relapse time between old and new relapse time
           ltrRL=max(tRLj,tRLjnew); % index of column for later relapse time between old and new relapse time
           erlrRLR=min(tRLRj,tRLRjnew); % index of column for earlier relapse treatment time between old and new relapse treatment time
           ltrRLR=max(tRLRj,tRLRjnew); % index of column for earlier relapse treatment time between old and new relapse treatment time
           
           % Calculate new infection pressure
           idx=[erlrRL+1:ltrRL,erlrRLR+1:ltrRLR]; % index of columns to change
           lambdaHH_new(:,idx)=rateHH*hnew(:,idx)+lambdaHHPA(:,idx)+lambdaHHA(:,idx)+pold(3); % new infection pressure on each HH
           lambda_new(:,idx)=lambdaHH_new(ib,idx); % new infection pressure on each individual
           
           % Calculate new log-likelihood terms
           LL1new=L1(S,lambda_new);
           LL2new=L2(lambda_new,pold(7),tEm);
           LL4new=L4(lambda_new,pold(7),tAm);
           LLnew=LL1new+LL2new+LL4new;
           LLold=LL1old+LL2old+LL4old;
           
           % Calculate Metropolis-Hastings acceptance probability
           log_ap=LLnew-LLold;
           
           if log_ap > log(rand) % accept new values if acceptance probability > rand
               tRL(j)=tRLnew(j); tRLR(j)=tRLRnew(j); h(m,:)=hnew(m,:); lambda(:,idx)=lambda_new(:,idx); lambdaHH(:,idx)=lambdaHH_new(:,idx);
               LL1old=LL1new; LL2old=LL2new;
               LL4old=LL4new; % keep updated info
               acc_RLNO=acc_RLNO+1;
           else % else reject and keep old values
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
       j=RLO(i); % get index of relapse case missing just relapse treatment time
       tRLRnew(j)=tRL(j)+nbinrnd(r0,p0)+1; % draw new relapse treatment time by drawing KA onset-to-treatment time and adding it to relapse time
       
       if tRLRnew(j)<=min(min(min(tEM(j),tP(j)),tD(j))-1,tmax) % calculate log-likelihood if new relapse treatment time is before emigration/PKDL onset/death/end of study
           tRLRj=tRLR(j);
           tRLRjnew=tRLRnew(j);
           m=(IPNIA==j); % index of case in infectiousness matrix
           hnew(m,tRLRj+1:tRLRjnew)=1; % add infectiousness if new relapse treatment time is later
           hnew(m,tRLRjnew+1:tRLRj)=0; % remove infectiousness if new relapse treatment time is earlier
           erlrRLR=min(tRLRj,tRLRjnew); % index of column for earlier relapse treatment time between old and new relapse treatment time
           ltrRLR=max(tRLRj,tRLRjnew); % index of column for later relapse treatment time between old and new relapse treatment time
           
           % Calculate new infection pressure
           idx=erlrRLR+1:ltrRLR; % index of columns to change
           lambdaHH_new(:,idx)=rateHH*hnew(:,idx)+lambdaHHPA(:,idx)+lambdaHHA(:,idx)+pold(3); % new infection pressure on each HH
           lambda_new(:,idx)=lambdaHH_new(ib,idx); % new infection pressure on each individual
           
           % Calculate new log-likelihood terms
           LL1new=L1(S,lambda_new);
           LL2new=L2(lambda_new,pold(7),tEm);
           LL4new=L4(lambda_new,pold(7),tAm);
           LLnew=LL1new+LL2new+LL4new;
           LLold=LL1old+LL2old+LL4old;
           
           % Calculate Metropolis-Hastings acceptance probability
           log_ap=LLnew-LLold;
           
           if log_ap > log(rand) % accept new values if acceptance probability > rand
               tRLR(j)=tRLRnew(j); h(m,:)=hnew(m,:); lambda(:,idx)=lambda_new(:,idx); lambdaHH(:,idx)=lambdaHH_new(:,idx);
               LL1old=LL1new; LL2old=LL2new;
               LL4old=LL4new; % keep updated info
               acc_RLO=acc_RLO+1;
           else % else reject and keep old values
               tRLRnew(j)=tRLR(j); hnew(m,:)=h(m,:); lambda_new(:,idx)=lambda(:,idx); lambdaHH_new(:,idx)=lambdaHH(:,idx); % keep old values, don't change log-likelihood
               rej_RLO=rej_RLO+1;
           end
       else % otherwise reject immediately
           tRLRnew(j)=tRLR(j);
           rej_RLO=rej_RLO+1;
       end
    end
    
    %% SAVE PARAMETER VALUES AND PLOT PROGRESS
    % Save parameters, log-likelihood, and presymptomatic infection times,
    % asymptomatic infection and recovery times and other missing data
    p(k+1,:)=pold;
    p1(k)=p1new;
    K0(k)=K0old;
    terms(k,:)=[LL1old,LL2old,LL3old,LL4old,LL5old,LL6old];
%     terms(k,:)=[LL1old,LL2old,LL3old,LL4old,LL5old];
    LL(k)=sum(terms(k,:));
    IPs(:,k)=IPold(I);
    tEs(:,k)=tE(I);
    tAs([SusA;PA],k)=tA([SusA;PA]);
    tRAs([SusA;PA],k)=tRA([SusA;PA]);
    tIsNONR(:,k)=tI(NONR);
    tRsNONR(:,k)=tR(NONR);
    tIsRNO(:,k)=tI(RNO);
    tRsONR(:,k)=tR(ONR);
    tIsANONR(:,k)=tI(ANONR);
    tRsANONR(:,k)=tR(ANONR);
    tRsAONR(:,k)=tR(AONR);
    tRLRsRLO(:,k)=tRLR(RLO);
    tRLsRLNO(:,k)=tRL(RLNO);
    tRLRsRLNO(:,k)=tRLR(RLNO);
    
    % Update empirical mean and covariance for proposal distribution
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
        
        % Plot output
        if plotOutpt && k>burnin % ignore burn-in period
            z=burnin+1:k; % iterations to plot
            figure(4);
            [mode_p,HPDI,mode_p1,HPDI1]=PlotOutput(z,LL,p,np,pname,priorp,p1,a,b,n,tmax,I,RpreD,DpreR,RL,tI,tR,tD,tRL,tRLR,nbins,scrnsz);            
            figure(5);
            PlotTrace(z,p,np,pname,p1,mode_p,HPDI,mode_p1,HPDI1,scrnsz)
            drawnow
            for j=1:np
                fprintf(['mode ' pname{j} '=%6.4g\n'], mode_p(j));
            end
            fprintf('mode p1=%6.4g\n', mode_p1);
        end
    end
end
 
%% DISPLAY ACCEPTANCE RATES
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

%% SAVE OUTPUT
% Remove first row (initial values) of p
p=p(2:end,:);

% Clear large variables
clear data dHH dHHsqrd KHH KHH_new rateHHA rateHHA_new rateHH rateHH_new

% Save workspace
save(rslts,'-v7.3')
