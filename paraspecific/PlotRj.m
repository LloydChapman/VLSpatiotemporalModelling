function PlotRj(Rj,idx,xlbl,ylbl)
figure;
meanRj=mean(Rj,1,'omitnan');
errorbar(idx,meanRj,meanRj-quantile(Rj,0.025,1),quantile(Rj,0.975,1)-meanRj,'.','MarkerSize',14)
set(gca,'FontSize',14)
xlabel(xlbl); ylabel(ylbl)
xlim([0 size(Rj,2)])