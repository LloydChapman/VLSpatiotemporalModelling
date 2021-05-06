function PlotAsxInfctnTimePstrDistn(tA,probA,tmax,idx)
histogram(tA,'Norm','prob','BinM','int'); hold on
plot(0:tmax+1,probA,'r','LineWidth',1)
set(gca,'FontSize',16)
xlabel(['$$A_{' num2str(idx) '}$$ (month)'],'Fontsize',18,'Interpreter','latex')
ylabel('Probability')
legend({['$$\mathrm{P}(A_{' num2str(idx) '}|\mathbf{X},\mathbf{Y},\theta)$$'],['$$q_{' num2str(idx) '}(A_{' num2str(idx) '})$$']},'Interpreter','latex','FontSize',18)
