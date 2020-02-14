function [M,h]=PlotInfctnTimePstrDistn(tE,tI,r,p,idx)
h=histogram(tE(~isnan(tE)),'Normalization','probability','BinMethod','integers'); hold on
% Find mode
BinCentres=(h.BinEdges(1:end-1)+h.BinEdges(2:end))/2;
[y,i]=max(h.Values);
M=BinCentres(i);
% Plot mode and prior
prior=nbinpdf(tI-BinCentres-1,r,p);
y=max(y,max(prior));
plot([M M],[0 1.2*y],'m','LineWidth',1.5)
% Plot prior distn
plot(BinCentres,prior,'g','LineWidth',1.5);
axis([min(h.BinEdges) max(h.BinEdges) 0 1.05*y])
set(gca,'FontSize',16);
xlabel(['$$E_{' num2str(idx) '}$$'],'Fontsize',16,'Interpreter','latex')
ylabel('Probability','FontSize',16)
