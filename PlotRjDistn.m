function PlotRjDistn(Rj,dim,xlbl,intrprtr,clr,varargin)
if nargin==6
    edges=varargin{1};
    hf=histogram(mean(Rj,dim,'omitnan'),edges);
else
    hf=histogram(mean(Rj,dim,'omitnan'));
end
hf.FaceColor=clr;
set(gca,'FontSize',14)
xlabel(xlbl,'Interpreter',intrprtr); ylabel('Count')