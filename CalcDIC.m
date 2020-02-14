function DIC=CalcDIC(rslts1,burnin1,ss,str)

rng=[];

load(rslts1)

if nargin==1
    burnin1=round(niters/10);
    ss='mode';
    str='DIC2';
end

% Overwrite asymptomatic infection and recovery times for PKDL cases w/o
% prior VL, if they were not saved, with values from final iteration
if tAs(PA(1))==tmax+2
    tAs(PA,:)=repmat(tA(PA),1,niters);
    tRAs(PA,:)=repmat(tRA(PA),1,niters);
end

% Load data
% load(db)
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

if ~exist('tLm','var')
    tLm=[];
end

[D,LLM]=dev(p,z,p1,nIPNIA,nIandP,IandP,I,tIM,tP,tRP,tEM,tmax,hv,nPI,PI,nI,tR,nIMP,IMP,nA,nIMI,I1,tI,NONR,tIsNONR,tRorD,tRsNONR,RNO,tIsRNO,ONR,tRsONR,ANONR,tIsANONR,tRsANONR,AONR,tRsAONR,RLO,tRLsRLO,tRLRsRLO,RLNO,tRLsRLNO,tRLRsRLNO,tEs,tAs,tRAs,IPs,n,IM_IN,IM_OUT,preB,preIM,tPm,tDm,tEMm,prevK,maxIP,IpreEXTIM,EXTIMsoonI,IpreINTIM,PpreINTIM,PpreEXTIM,h0,A,IMI,nRL,RL,tRL,tRLR,tD,INTMIG_OUT,rng,nPA,PA,dHH,dHHsqrd,typ,nHH,ib,f,u,d0,IPNIA,r1,age,S0,p2,S0PA,actvAPA,prevAPA,ss);
MD=-2*mean(LL(z));
DIC=dic(D,MD);

save([str '_' rslts1],'D','LLM','MD','DIC')