function hf=PlotMCTAndCI(x,MCT,CI,clr,alph,xlbl,ylbl)
x2=[x,fliplr(x)];
hf=plot(x,MCT,'LineWidth',1,'Color',clr); hold on
CIR=[CI(1,:),fliplr(CI(2,:))];
hf1=fill(x2,CIR,clr,'LineStyle','none');
hf1.FaceAlpha=alph;
set(gca,'FontSize',14)
xlabel(xlbl)
ylabel(ylbl,'Interpreter','latex')
xlim([x(1) x(end)])
hold off