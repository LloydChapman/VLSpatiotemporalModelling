function [meanmeandists,meanmeantimes]=PlotDistsAndTimes(infctr,src,dists,times,i,hf1,hf2,clr)
nsmpls=size(infctr,2);
infctri=infctr(src==i);
unqinfctr=unique(infctri);
meandists=NaN(numel(unqinfctr),nsmpls);
meantimes=NaN(numel(unqinfctr),nsmpls);
for j=1:nsmpls   
    infctrj=infctr(src(:,j)==i,j);
    distsj=dists(src(:,j)==i,j);
    timesj=times(src(:,j)==i,j);
    for k=1:numel(unqinfctr)
        if any(infctrj==unqinfctr(k))
            meandists(k,j)=mean(distsj(infctrj==unqinfctr(k)));
            meantimes(k,j)=mean(timesj(infctrj==unqinfctr(k)));
        else
            meandists(k,j)=NaN;
            meantimes(k,j)=NaN;
        end
    end   
end
meanmeandists=mean(meandists,2,'omitnan');
meanmeantimes=mean(meantimes,2,'omitnan');
edges1=0:50:max(meanmeandists)+50;
edges2=1:max(meanmeantimes)+1;
figure(hf1); histogram(meanmeandists,edges1,'FaceColor',clr)
figure(hf2); histogram(meanmeantimes,edges2,'FaceColor',clr) %'BinLimits',[1 max(meanmeantimes)],

