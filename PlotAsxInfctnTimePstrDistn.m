function PlotAsxInfctnTimePstrDistn(tA,probA,tmax,idx)
histogram(tA,'Norm','prob','BinM','int'); hold on
plot(0:tmax+1,probA,'r','LineWidth',1)
set(gca,'FontSize',14)
xlabel(['$$A_{' num2str(idx) '}$$ (months)'],'Fontsize',16,'Interpreter','latex')
ylabel('Probability')
legend({['$$\mathrm{P}(A_{' num2str(idx) '}|\mathbf{X},\mathbf{Y},\theta)$$'],['$$q_{' num2str(idx) '}(A_{' num2str(idx) '})$$']},'Interpreter','latex')

% h=histogram(tE(~isnan(tE)),'Normalization','probability','BinMethod','integers'); hold on
% % Find mode
% BinCentres=(h.BinEdges(1:end-1)+h.BinEdges(2:end))/2;
% [y,i]=max(h.Values);
% M=BinCentres(i);
% % Plot mode and prior
% prior=nbinpdf(tI-BinCentres-1,r,p);
% y=max(y,max(prior));
% plot([M M],[0 1.2*y],'m','LineWidth',1.5)
% % Plot prior distn
% plot(BinCentres,prior,'g','LineWidth',1.5);
% axis([min(h.BinEdges) max(h.BinEdges) 0 1.05*y])
% set(gca,'FontSize',16);
% xlabel(['$$E_{' num2str(idx) '}$$'],'Fontsize',16,'Interpreter','latex')
% ylabel('Density','FontSize',16)
