function [infctr,src,infctrmax,srcmax,RjA,RjI,RjP,Rj,Rts,Rt,diA,diI,diP,di,tiA,tiI,tiP,ti]=GetInfctrAndSrc(infctn,rateA,rate,ratePA,hA,h,hP,hPA,p,I1,A1,nI1,tA,onsetI,onsetP,onsetPA,dHH,ib,PA,IPNIA,maxIP,tmax,stat,n,nIPNIA)
%% EFFECTIVE REPRODUCTION NUMBERS
rateI1=rate(I1,:);
lambdaAij=NaN(nI1,n);
lambdaAij(:,A1)=rateA(I1,A1).*hA(A1,infctn)';
lambdaIij=rateI1.*h(:,infctn)';
lambdaPij=rateI1.*hP(:,infctn)';
lambdaPAij=ratePA(I1,:).*hPA(:,infctn)';
lambdaij=[lambdaAij,lambdaIij,lambdaPij,lambdaPAij,p(3)*ones(nI1,1)];
% lambdai=sum(lambdaij,2);
% pi=sum(1-exp(-lambdaij),2,'omitnan');
lambdai=sum(lambdaij,2,'omitnan');
% onset=[tA;onsetI;onsetP;onsetPA;NaN];
onset=[tA;onsetI;onsetPA;NaN];

% Rts=zeros(3,tmax);
dAI1=dHH(ib(I1),ib);
dI1=dHH(ib(I1),ib(IPNIA));
dPAI1=dHH(ib(I1),ib(PA));
tA(tA==tmax+1)=NaN;
tijA=bsxfun(@minus,infctn,double(tA)');
tijI=bsxfun(@minus,infctn,onsetI');
tijP=bsxfun(@minus,infctn,onsetP');
tijPA=bsxfun(@minus,infctn,onsetPA');
[RjA,diA,tiA]=CalcRjAnddi(lambdaAij,lambdai,dAI1,tijA);
% Rts(1,:)=CalcRt(RjA,onsetA,tmax);
[RjI,diI,tiI]=CalcRjAnddi(lambdaIij,lambdai,dI1,tijI);
% Rts(2,:)=CalcRt(RjI,onsetI,tmax);
[RjP,diP,tiP]=CalcRjAnddi([lambdaPij,lambdaPAij],lambdai,[dI1,dPAI1],[tijP,tijPA]);
% Rts(3,:)=CalcRt(RjP,[onsetP;onsetPA],tmax);
% RjPA=CalcRj(lambdaPAij,lambdai);
% RtPA=CalcRt(RjPA,onsetPA,tmax);
[Rj,di,ti,pij_ti]=CalcRjAnddi(lambdaij,lambdai,[dAI1,dI1,dI1,dPAI1,NaN(nI1,1)],[tijA,tijI,tijP,tijPA,NaN(nI1,1)]);
% Add together numbers of secondary infections from VL and PKDL episodes of 
% cases
Rj=[Rj(:,1:n),Rj(:,n+1:n+nIPNIA)+Rj(:,n+nIPNIA+1:n+2*nIPNIA),Rj(:,n+2*nIPNIA+1:size(Rj,2))];
[Rt,Rts]=CalcRt(Rj,RjA,RjI,RjP,onset,tA,onsetI,maxIP,tmax);

%% INFECTORS AND INFECTION SOURCES
% Make index vector of infectors and matrix of their statuses at
% infection times
IPNIAbckgrnd=[(1:n)';IPNIA;IPNIA;PA;n+1];
statbckgrnd=[stat(:,infctn)',9*ones(nI1,1)];

%% Infector based on posterior distribution
% Draw infector from posterior distn of possible infectors
infctr=IPNIAbckgrnd(sum(rand(nI1,1)>=cumsum(pij_ti,2,'omitnan'),2)+1);
src=statbckgrnd(sub2ind([nI1 n+1],(1:nI1)',infctr));

%% Most likely infector
% Find maximum infectious pressures (incl. background) on KA cases at
% their infection times
[~,i]=max(lambdaij,[],2);
% Find and store index of most likely infector
infctrmax=IPNIAbckgrnd(i);

% Find most likely infection source (asymptomatic, pre-symptomatic, KA,
% PKDL, background)
srcmax=statbckgrnd(sub2ind([nI1 n+1],(1:nI1)',infctrmax));