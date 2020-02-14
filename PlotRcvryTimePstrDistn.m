function [M,h]=PlotRcvryTimePstrDistn(tR,tI,r,p,idx)
h=histogram(tR,'Normalization','probability','BinMethod','integers'); hold on
% Find mode
BinCentres=(h.BinEdges(1:end-1)+h.BinEdges(2:end))/2;
[y,i]=max(h.Values);
M=BinCentres(i);
% Plot mode and prior
prior=nbinpdf(BinCentres-tI-1,r,p);
y=max(y,max(prior));
plot([M M],[0 10*y],'m','LineWidth',1.5)
% Plot prior distn
plot(BinCentres,prior,'g','LineWidth',1.5);
axis([min(h.BinEdges) max(h.BinEdges) 0 y])
set(gca,'FontSize',16);
xlabel(['$$R_{' num2str(idx) '}$$'],'Fontsize',13,'Interpreter','latex')
ylabel('Density')
