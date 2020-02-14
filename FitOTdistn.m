function [r,p]=FitOTdistn(tI,tR)
% Calculate onset-to-treatment (OT) intervals 
OT=tR-tI;
% Remove missing OT times
OT=OT(~isnan(OT));
% Fit NB distn to OT-1 times (subtract 1 as support for NB starts from 0 in MATLAB)
pars=nbinfit(OT-1);
r=pars(1);
p=pars(2);
% % Plot fit
% figure; histogram(OT,'Normalization','pdf','BinMethod','integers'); hold on
% x=1:max(OT); plot(x,nbinpdf(x-1,r,p),'LineWidth',1); xlabel('OT (mnths)'); ylabel('Probability'); hold off
% legend('data','g(OT)');
% % saveas(gcf,'BNVL_OTdistn');
% % saveas(gcf,'BNVL_OTdistn.eps','epsc');