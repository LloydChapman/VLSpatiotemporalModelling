function LL=RecalclogL(n,tmax,tE,I2,I1,tA,tRA,IM_IN,IM_OUT,preB,preIM,tDm,tEMm,prevK,maxIP,IpreEXTIM,EXTIMsoonI,IpreINTIM,PpreINTIM,PpreEXTIM,h,nI,I,tIM,tI,h0,tEM,tRorD,nA,nPI,A,nIMI,IMI,nRL,RL,IPNIA,tRL,tRLR,tD,h4,INTMIG_OUT,rng,nPA,PA,tP,tRP,hv,dHH,dHHsqrd,alpha,beta,typ,nHH,ib,f,u,delta,d0,epsilon,pI,IP,r1,p1,age,S0,lambda0,p2,S0PA,actvAPA,prevAPA)
    
tEm=false(n,tmax);
tEm((tE(I2)-1)*n+uint32(I1))=1;
    
Asx=uint32(find(tA>=1 & tA<=tmax));
tAm=false(n,tmax);
tAm((tA(Asx)-1)*n+Asx)=1;
prevA=uint32(find(tA==0 & tRA==0));
actvA=uint32(find(tA==0 & tRA>0));
IM_INprevAactvA=IM_IN(ismember(IM_OUT,[prevA;actvA]));

S=1-max(preB,preIM)-max(max(max(cumsum(tEm,2),cumsum(tAm,2)),cumsum(tDm,2)),cumsum(tEMm,2)); % don't remove LST+ individuals
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

for i=1:nI
    j=I(i);
    h(i,max(max(0,tIM(j)),tE(i))+1:tI(j))=h0;
    h(i,max(tIM(j),tI(j))+1:min(tEM(j),tRorD(j)))=1;
end
for i=1:nA
    h(nI+nPI+i,1:tRorD(A(i)))=1;
end
for i=1:nIMI
    j=IMI(i);
    h(nI+nPI+nA+i,tIM(j)+1:tRorD(j))=1;
end
for i=1:nRL
    j=RL(i);
    h(IPNIA==j,tRL(j)+1:min(min(tRLR(j),tEM(j)),tmax))=1;
end
hA=zeros(n,tmax);
A1=find(tA>=0 & tA<tmax+1 & tRA>0 & isnan(tP));
for i=1:numel(A1)
    j=A1(i);
    hA(j,tA(j)+1:min(min(min(tRA(j),tEM(j)),tD(j)),tmax))=h4;
end
RAobs1actvA=find(tA==0 & tRA>=tEM & INTMIG_OUT);
RAobs2actvA=IM_IN(ismember(IM_OUT,RAobs1actvA));
RAobs1=find(tA>0 & tA<tEM & tRA>=tEM & INTMIG_OUT);
RAobs2=IM_IN(ismember(IM_OUT,RAobs1));
for i=1:numel(RAobs2actvA)
    j=RAobs1actvA(i);
    j1=RAobs2actvA(i);
    hA(j1,rng(j,2)+1:min(min(tRA(j),rng(j1,2)),tmax))=h4;
end
for i=1:numel(RAobs2)
    j=RAobs1(i);
    j1=RAobs2(i);
    hA(j1,rng(j,2)+1:min(min(tRA(j),rng(j1,2)),tmax))=h4;
end
hA=sparse(hA);
hPA=zeros(nPA,tmax);
for i=1:nPA
    j=PA(i);
    hPA(i,max(1,tA(j)+1):tRA(j))=h4;
    hPA(i,tP(j)+1:min(min(tRP(j),tEM(j)),tmax))=hv(j);
end

KHH=Knl_fast(dHH,dHHsqrd,alpha,beta,typ,n,nHH,ib,f,1:n);
rateHHA=beta*KHH;
if ismember(4,u)
    rateHHA=rateHHA+delta*d0;
end
rateHH=rateHHA(:,ib(IPNIA));
rateHHPA=rateHHA(:,ib(PA));
%     lambdaHHA=rateHHA(:,[A1,RAobs2actvA,RAobs])*hA([A1,RAobs2actvA,RAobs],:);
%     lambdaHHPA=rateHHPA*hPA;
lambdaHH=rateHH*h+rateHHA(:,ib([A1;RAobs2actvA;RAobs2]))*hA([A1;RAobs2actvA;RAobs2],:)+rateHHPA*hPA+epsilon;
lambda=lambdaHH(ib,:);

LL1=L1(S,lambda);
LL2=L2(lambda,pI,tEm);
LL3=L3(IP,r1,p1);
LL4=L4(lambda,pI,tAm);
LL5=L5(age,S0,actvA,prevA,lambda0,pI,p2);
LL6=L5(age,S0PA,actvAPA,prevAPA,lambda0,pI,p2);
LL=LL1+LL2+LL3+LL4+LL5+LL6;
