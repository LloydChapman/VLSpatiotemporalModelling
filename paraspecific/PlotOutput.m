function [mode_p,HPDI,mode_p1,HPDI1]=PlotOutput(z,LL,p,np,pname,priorp,p1,a,b,n,tmax,I,RpreD,DpreR,RL,tI,tR,tD,tRL,tRLR,nbins,scrnsz)
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
tIm=false(n,tmax);
tIm((tI(I)-1)*n+I)=1;
tRorDm=false(n,tmax);
tRorDm((tR(RpreD)-1)*n+RpreD)=1;
tRorDm((tD(DpreR)-1)*n+DpreR)=1;
tRLm=false(n,tmax);
tRLm((tRL(RL)-1)*n+RL)=1;
tRLRm=false(n,tmax);
tRLRm((tRLR(RL)-1)*n+RL)=1;
NI=sum(cumsum(tIm,2)-cumsum(tRorDm,2)+cumsum(tRLm,2)-cumsum(tRLRm,2));
plot(NI,'r');
xlim([1 Inf])
xlabel('t (months)');
ylabel('No. of VL cases');