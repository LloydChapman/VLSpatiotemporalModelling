function PlotRjDistn(Rj,dim,xlbl,intrprtr,clr,varargin)
if nargin==6
    edges=varargin{1};
    hf=histogram(mean(Rj,dim,'omitnan'),edges,'Normalization','pdf');
else
    hf=histogram(mean(Rj,dim,'omitnan'),'Normalization','pdf');
end
hf.FaceColor=clr;
set(gca,'FontSize',14)
xlabel(xlbl,'Interpreter',intrprtr); ylabel('Density')