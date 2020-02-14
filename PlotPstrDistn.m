function [M,HPDI]=PlotPstrDistn(p,pname,nbins,varargin)
h=histogram(p,nbins,'Normalization','pdf'); hold on
y=max(h.Values);
if all(p==p(1))
    M=p(1);
    HPDI=[p(1) p(1)];
else
    if nargin==6
        binmthd=varargin{3};
        [M,HPDI]=CalcModeAndHPDI(p,nbins,binmthd);
    else
        [M,HPDI]=CalcModeAndHPDI(p,nbins);
    end
end
% Plot HPDI
plot([HPDI(1) HPDI(1)],[0 1.05*y],'r',[HPDI(2) HPDI(2)],[0 1.05*y],'r','LineWidth',1.5)
% indcs=min(indcs):max(indcs); % make sure range for HPDI is continuous
% bar((hp.BinEdges(indcs)+hp.BinEdges(indcs+1))/2,hp.Values(indcs),1,'FaceColor',[0.85 0.325 0.098])
% Plot mode
plot([M M],[0 1.05*y],'m','LineWidth',1.5)
% Plot prior distn
if nargin>=5
    prior=varargin{1};
    priorp=varargin{2};
    if numel(priorp)==1
        plot(h.BinEdges,pdf(prior,h.BinEdges,priorp),'g','LineWidth',1.5)
    elseif numel(priorp)==2
        plot(h.BinEdges,pdf(prior,h.BinEdges,priorp(1),priorp(2)),'g','LineWidth',1.5)
    elseif numel(priorp)==3
        plot(h.BinEdges,pdf(prior,h.BinEdges,priorp(1),priorp(2),priorp(3)),'g','LineWidth',1.5)
    end
end
axis([min(h.BinEdges) max(h.BinEdges) 0 1.05*y])
% set(gca,'Fontsize',16)
if ismember(pname,{'beta','alpha','epsilon','delta','lambda_0'})
    xlabel(['\' pname],'FontSize',13)
%     xlabel(['\' pname],'FontSize',20)
else
    xlabel(['$$' pname '$$'],'FontSize',13,'Interpreter','latex')
%     xlabel(pname,'FontSize',16);
end
ylabel('Density')
% ylabel('Probability density','FontSize',14)
hold off
% saveas(gcf,['PosteriorDistn' pname{i} '.eps'],'epsc')