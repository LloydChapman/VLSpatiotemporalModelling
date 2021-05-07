function PlotTrace(z,p,np,pname,p1,mode_p,HPDI,mode_p1,HPDI1,scrnsz)
set(gcf, 'Position', [round(scrnsz(3)/2) 70 round(scrnsz(3)/2) scrnsz(4) - 150]);
z1=[z(1) z(end)];
for j=1:np
    subplot(floor(np/2)+1, 2, j)
    plot(z,p(z, j)); hold on
    plot(z1,[mode_p(j) mode_p(j)],'-m',z1,[HPDI(j,1) HPDI(j,1)],'r-',z1,[HPDI(j,2) HPDI(j,2)],'r-');
    axis([z1 -Inf Inf]);
    xlabel('Iteration')
    if ismember(pname{j},{'beta','alpha','epsilon','delta','lambda_0'})
        ylabel(['\' pname{j}],'FontSize',13,'Rot',0)       
    else
        ylabel(['$$' pname{j} '$$'],'FontSize',13,'Rot',0,'Interpreter','latex')
    end
    hold off
end
subplot(floor(np/2)+1,2,np+1)
plot(z,p1(z)); hold on
plot(z1,[mode_p1 mode_p1],'-m',z1,[HPDI1(1) HPDI1(1)],'r-',z1,[HPDI1(2) HPDI1(2)],'r-');
axis([min(z) max(z) -Inf Inf]);
xlabel('Iteration')
ylabel('$$p$$','Fontsize',13,'Rot',0,'Interpreter','latex')
hold off