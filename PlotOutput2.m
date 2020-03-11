function [mode_p,HPDI,mode_p1,HPDI1]=PlotOutput2(z,LL,p,np,pname,priorp,p1,a,b,n,tmax,I,RpreD,DpreR,OR,NONR,RNO,ONR,A,ANONR,AONR,RL,RLO,RLNO,tI,tR,tD,tRL,tRLR,tIsNONR,tIsRNO,tRsNONR,tRsONR,tIsANONR,tRsANONR,tRsAONR,tRLsRLNO,tRLRsRLO,tRLRsRLNO,nbins,scrnsz)
set(gcf, 'Position', [0 70 round(scrnsz(3)/2) scrnsz(4) - 150]);
mode_p=zeros(1,np); % modal parameter values
HPDI=zeros(np,2); % highest posterior density intervals
subplot(3+floor(np/2), 2, [1 2])
plot(z,LL(z));
axis([z(1) z(end) min(LL(z)) max(LL(z))]);
xlabel('Iteration');
ylabel('Log likelihood');
for j=1:np    
    subplot(3+floor(np/2), 2, 2+j)
    [mode_p(j),HPDI(j,:)]=PlotPstrDistn(p(z,j),pname{j},nbins,'gamma',priorp{j});
end
subplot(3+floor(np/2),2,np+3)
[mode_p1,HPDI1]=PlotPstrDistn(p1(z),'p',nbins,'beta',[a,b]);
subplot(3+floor(np/2),2,[5+2*floor(np/2) 6+2*floor(np/2)])

% Number of iterations to plot
npp=min(100,numel(z));
pp=randi([z(1) z(end)],1,npp);

% Calculate numbers of KA cases over time for different MCMC iterations
tIs=repmat(tI,1,npp);
tIs(NONR,:)=tIsNONR(:,pp);
tIs(RNO,:)=tIsRNO(:,pp);
tIs(ANONR,:)=tIsANONR(:,pp);
tIm=false(n,tmax,npp);
tRs=repmat(tR,1,npp);
tRs(NONR,:)=tRsNONR(:,pp);
tRs(ONR,:)=tRsONR(:,pp);
tRs(ANONR,:)=tRsANONR(:,pp);
tRs(AONR,:)=tRsAONR(:,pp);
tRorDm=false(n,tmax,npp);
tRLs=repmat(tRL,1,npp);
tRLs(RLNO,:)=tRLsRLNO(:,pp);
tRLm=false(n,tmax,npp);
tRLRs=repmat(tRLR,1,npp);
tRLRs(RLO,:)=tRLRsRLO(:,pp);
tRLRs(RLNO,:)=tRLRsRLNO(:,pp);
tRLRm=false(n,tmax,npp);
for i=1:npp
tIm((i-1)*n*tmax+(tIs(I,i)-1)*n+I)=1;
tIm(A(tRs(A,i)>=1),1,i)=1;
tRorDm((i-1)*n*tmax+(tRs(RpreD,i)-1)*n+RpreD)=1;
tRorDm((i-1)*n*tmax+(tD(DpreR)-1)*n+DpreR)=1;
tRorDm((i-1)*n*tmax+(tRs(A(tRs(A,i)>=1),i)-1)*n+A(tRs(A,i)>=1))=1;
tRLm((i-1)*n*tmax+(tRLs(RL,i)-1)*n+RL)=1;
tRLRm((i-1)*n*tmax+(tRLRs(RL,i)-1)*n+RL)=1;
end
numI=sum(cumsum(tIm,2)-cumsum(tRorDm,2)+cumsum(tRLm,2)-cumsum(tRLRm,2));

% Calculate numbers of KA cases over time with observed onset and recovery
tIm1=false(n,tmax);
tIm1((tI(OR)-1)*n+OR)=1;
ORpreD=setdiff(OR,DpreR);
ODpreR=intersect(OR,DpreR);
tRorDm1=false(n,tmax);
tRorDm1((tR(ORpreD)-1)*n+ORpreD)=1;
tRorDm1((tD(ODpreR)-1)*n+ODpreR)=1;
numI1=sum(cumsum(tIm1,2)-cumsum(tRorDm1,2)+cumsum(tRLm,2)-cumsum(tRLRm,2));

% Plot imputed case numbers vs observed case numbers
jj=1:tmax;
plot(jj,squeeze(numI(1,jj,:)),'r'); hold on
plot(jj,numI1(jj),'k--','LineWidth',2); hold off
xlim([jj(1) Inf])
xlabel('t (months)');
ylabel('No. of VL cases');