clear

%% Load data
load('C:\Users\timpo\OneDrive - University of Warwick\taubern_baybern\raw_data_plus_cleaning\matlab_bayesianmodel\data_final2.mat');
data.MOS_RX_NEW_SX = str2double(data.MOS_RX_NEW_SX);
% Number of individuals
n=size(data,1);

% Duration of infectiousness after starting treatment - assume cases stop being infectious shortly after commencing treatment
durRX=0; % months

% Set start year and month
startyr=2002;
startmo=1;
 
% Set end year and month
endyr=2010;
endmo=12;
 
% Set maximum incubation period for cases with onset at start of study
maxIP=12; %6; %11; %

% Set maximum incubation period for cases with onset shortly after
% migration in
maxIP_IM=4;

%% Make vectors of event times
% Make STATA month origin
origin=stata_month(startyr,startmo)-1;

% Make indicator variables for 1st and 2nd observations for individuals who internally migrated 
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
tR=data.KARX+durRX-origin;
% Add 1 month of infectiousness for treatment failure KA cases
tR(RXF)=tR(RXF)+1;
tRL=data.KARX+data.MOS_RX_NEW_SX-origin;
tRL(RXF)=NaN;
tRLR=NaN(n,1);
tP=data.PKDL-origin;
tRP=data.PKDL+data.PKDL_DUR-origin;
% Overwrite resolution times with treatment times for treated PKDL cases.
tRP(RXP)=max(data.PKRX(RXP),data.PKRX2(RXP))-origin; % use 2nd treatment time for case with 2 PKDL treatments
tB=data.DOB-origin;
tD=data.DEATH-origin;
tIM=data.MIG_IN-origin;
tEM=data.MIG_OUT-origin;

RLO=find(REL&~isnan(tRL));

%% Calculate fixed parameters for simulations
% Fit negative binomial distribution to onset-to-treatment (OT) times
[r0,p0]=FitOTdistn(tI(~KothrObs),tR(~KothrObs));
OT=tR(~KothrObs)-tI(~KothrObs);
% Plot fit
figure; histogram(OT(~isnan(OT)),'Normalization','probability','BinMethod','integers'); hold on
x=1:max(OT); plot(x,nbinpdf(x-1,r0,p0),'LineWidth',1)
set(gca,'FontSize',14)
xlabel('VL onset-to-treatment time (months)'); ylabel('Probability')
xlim([0.5 Inf])
legend('data','fitted pmf')
hold off
saveas(gcf,'VL_OTdistn')
saveas(gcf,'VL_OTdistn.eps','epsc') 

% Fit negative binomial distribution to KA-treatment-to-PKDL-onset times
RP=tP-tR;
pars1=nbinfit(RP(RP>=0 & ~PothrObs));
r3=pars1(1);
p3=pars1(2);
% Plot fit
figure; histogram(RP(RP>=0 & ~PothrObs),'Normalization','probability','BinMethod','integers'); hold on
x=0:max(RP); plot(x,nbinpdf(x,r3,p3),'LineWidth',1)
set(gca,'FontSize',14)
xlim([-0.5 Inf])
xlabel('VL-treatment-to-PKDL-onset time (months)'); ylabel('Probability')
legend('data','fitted pmf')
hold off
saveas(gca,'VLRXtoPKDLdistn')
saveas(gca,'VLRXtoPKDLdistn.eps','epsc')

% Fit geometric distribution to KA-treatment-to-relapse times
p4=mle(tRL(RLO)-tR(RLO)-1,'distribution','geo');

% Fit negative binomial distribution to PKDL-onset-to-recovery times
pars2=nbinfit(tRP(~isnan(tRP)&~PothrObs)-tP(~isnan(tRP)&~PothrObs)-1);
r5=pars2(1);
p5=pars2(2);

% Save parameters
save('FixedParams','r0','p0','r3','p3','p4','r5','p5')





 