function [D,LL]=dev(p,z,p1,nIPNIA,nIandP,IandP,I,tIM,tP,tRP,tEM,tmax,hv,nPI,PI,nI,tR,nIMP,IMP,nA,nIMI,I1,tI,NONR,tIsNONR,tRorD,tRsNONR,RNO,tIsRNO,ONR,tRsONR,ANONR,tIsANONR,tRsANONR,AONR,tRsAONR,RLO,tRLRsRLO,RLNO,tRLsRLNO,tRLRsRLNO,tEs,tAs,tRAs,IPs,n,IM_IN,IM_OUT,preB,preIM,tDm,tEMm,prevK,maxIP,IpreEXTIM,EXTIMsoonI,IpreINTIM,PpreINTIM,PpreEXTIM,h0,A,IMI,nRL,RL,tRL,tRLR,tD,INTMIG_OUT,rng,nPA,PA,dHH,dHHsqrd,typ,nHH,ib,f,u,d0,IPNIA,r1,age,S0,p2,S0PA,actvAPA,prevAPA,ss)

if strcmp(ss,'mean')
    beta=mean(p(z,1));
    alpha=mean(p(z,2));
    epsilon=mean(p(z,3));
    delta=mean(p(z,4));
    lambda0=mean(p(z,5));
    h4=mean(p(z,6));
    pI=mean(p(z,7));
    p1=mean(p1(z));
elseif strcmp(ss,'mode')
    nbins=50;
    beta=CalcModeAndHPDI(p(z,1),nbins);
    alpha=CalcModeAndHPDI(p(z,2),nbins);
    epsilon=CalcModeAndHPDI(p(z,3),nbins);
    delta=CalcModeAndHPDI(p(z,4),nbins);
    lambda0=CalcModeAndHPDI(p(z,5),nbins);
    h4=CalcModeAndHPDI(p(z,6),nbins);
    pI=CalcModeAndHPDI(p(z,7),nbins);
    p1=CalcModeAndHPDI(p1(z),nbins);
end

h=zeros(nIPNIA,tmax);
for i=1:nIandP
    j=IandP(i);
    h(I==j,max(tIM(j),tP(j))+1:min(min(tRP(j),tEM(j)),tmax))=hv(j);
end
for i=1:nPI
    j=PI(i);
    h(nI+i,tP(j)+1:min(min(tRP(j),tEM(j)),tmax))=hv(j);
end
% Need to remove PKDL infectiousness for individual with simultaneous KA
% and PKDL here so we don't double count their contribution, as we treat 
% them as having KA infectiousness until they were treated for KA in MCMC
% code
PpreR=find(tR>tP&~isinf(tR));
h(ismember(I,PpreR),tP(PpreR)+1:tR(PpreR))=0;
for i=1:nIMP
    j=IMP(i);
    h(nI+nPI+nA+nIMI+i,tIM(j)+1:min(tRP(j),tmax))=hv(j);
end
I2=ismember(I,I1);

LL=zeros(numel(z),1);
for i=1:numel(z)
%     i
    k=z(i);
    
    tI(NONR)=tIsNONR(:,k);
    tRorD(NONR)=tRsNONR(:,k);
    tI(RNO)=tIsRNO(:,k);
    tRorD(ONR)=tRsONR(:,k);
    tI(ANONR)=tIsANONR(:,k);
    tRorD(ANONR)=tRsANONR(:,k);
    tRorD(AONR)=tRsAONR(:,k);
    tRLR(RLO)=tRLRsRLO(:,k);
    tRL(RLNO)=tRLsRLNO(:,k);
    tRLR(RLNO)=tRLRsRLNO(:,k);
    
    tE=uint32(tEs(:,k));
    tA=uint32(tAs(:,k));
    tRA=uint32(tRAs(:,k));
    IP=IPs(:,k);
    
    LL(i)=RecalclogL(n,tmax,tE,I2,I1,tA,tRA,IM_IN,IM_OUT,preB,preIM,tDm,tEMm,prevK,maxIP,IpreEXTIM,EXTIMsoonI,IpreINTIM,PpreINTIM,PpreEXTIM,h,nI,I,tIM,tI,h0,tEM,tRorD,nA,nPI,A,nIMI,IMI,nRL,RL,IPNIA,tRL,tRLR,tD,h4,INTMIG_OUT,rng,nPA,PA,tP,tRP,hv,dHH,dHHsqrd,alpha,beta,typ,nHH,ib,f,u,delta,d0,epsilon,pI,IP,r1,p1,age,S0,lambda0,p2,S0PA,actvAPA,prevAPA);
end

D=-2*mean(LL);