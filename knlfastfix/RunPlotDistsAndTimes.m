function [meandists,meantimes]=RunPlotDistsAndTimes(infctr,src,dists,times,srcidx,clrs,name1,name2)
hf1=figure; hold on
hf2=figure; hold on
meandists=cell(numel(srcidx),1);
meantimes=cell(numel(srcidx),1);
for i=1:numel(srcidx)
    j=srcidx(i);
    [meandists{i},meantimes{i}]=PlotDistsAndTimes(infctr,src,dists,times,j,hf1,hf2,clrs(i,:));
end
ax1=get(hf1,'CurrentAxes');
set(ax1,'FontSize',14)
xlabel(ax1,'Mean distance of infectees from infector (m)'); ylabel(ax1,'Density')
legend(ax1,'VL','PKDL')
saveas(hf1,name1)
saveas(hf1,[name1 '.eps'],'epsc')
saveas(hf1,[name1 '.png'])
ax2=get(hf2,'CurrentAxes');
set(ax2,'FontSize',14)
xlabel(ax2,'Mean infector-onset-to-infectee-infection time (months)'); ylabel(ax2,'Density')
legend(ax2,'VL','PKDL')
saveas(hf2,name2)
saveas(hf2,[name2 '.eps'],'epsc')
saveas(hf2,[name2 '.png'])