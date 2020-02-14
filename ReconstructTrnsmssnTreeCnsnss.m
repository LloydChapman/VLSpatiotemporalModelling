% function ReconstructTrnsmssnTreeCnsnss(infctn,infctr,src,typinfctn,I1,dHH,ib,tI,tP,I,P,infctnA)
% function [gen,dist,time,chain]=ReconstructTrnsmssnTreeCnsnss(n,nI1,tmax,infctn,infctr,src,I1,dHH,ib,tI,tP,I,P,infctnA)
function [gen,dist,time,chain]=ReconstructTrnsmssnTreeCnsnss(rslts,infctn,infctr,src,infctnA)

db='';
load(rslts)
% load(db)
load('data_final2.mat')

% Select data for para
data=data(ismember(data.PARA,para),:);
% Rename longitude and latitude variables
data.Properties.VariableNames{'HHNEWLNG'}='longitude';
data.Properties.VariableNames{'HHNEWLAT'}='latitude';
% Calculate HH distance matrix
dHH=CalcHHDists(data);

nI1=numel(I1);
nsmpls=size(infctn,2);
gen=NaN(nI1+n,nsmpls);
dist=NaN(nI1+n,nsmpls);
time=NaN(nI1+n,nsmpls);
chain=cell(nI1+n,nsmpls);

for i=1:nsmpls
%     [typinfctnsrt(:,i),I1srt(:,i),infctnsrt(:,i),infctrsrt(:,i),srcsrt(:,i),gen(:,i),dist(:,i),time(:,i),chain(:,i)]=ReconstructTrnsmssnTree(infctn(:,i),infctr(:,i),src(:,i),typinfctn,I1,dHH,ib,tI,tP,I,P,infctnA(:,i));
    [gen(:,i),dist(:,i),time(:,i),chain(:,i)]=ReconstructTrnsmssnTree(n,nI1,tmax,infctn(:,i),infctr(:,i),src(:,i),I1,dHH,ib,tI,tP,I,P,infctnA(:,i));
end

figure;
boxplot(dist(:),gen(:))
set(gca,'FontSize',14)
xlabel('Infection generation')
ylabel('Distance from index case (m)')
% saveas(gcf,'DistFromIndexCaseVsGen.png')

figure;
boxplot(time(:),gen(:))
set(gca,'FontSize',14)
xlabel('Infection generation')
ylabel('Time since index case onset (months)')
% saveas(gcf,'TimeFromIndexCaseVsGen.png')

% cmap=distinguishable_colors(max(gen(:))+3);
% figure; axes('colororder',cmap(4:end,:)); hold on
% cmap=linspecer(max(gen(:)));
cmap=colormap(jet(max(gen(:))));
% cmap=cmap([5 1 9 13 2 10 4 3 6 7 11 8 12],:);
cmap=cmap(randperm(max(gen(:))),:);
figure; axes('colororder',cmap); hold on
for i=1:max(gen(:))
    lh(i)=plot(dist(gen==i),time(gen==i),'.');
    lt{i}=num2str(i);
end
set(gca,'FontSize',14)
xlabel('Distance from index case (m)')
ylabel('Time since index case onset (months)')
[lgd,icons]=legend(lh,lt);
% title(lgd,'Generation')
icons=findobj(icons,'Type','line');
icons=findobj(icons, '-not', 'Marker', 'none');
set(icons,'MarkerSize',16)
hold off
% saveas(gcf,'TimeVsDistFromIndexCase.png')

dist1=dist; dist1(mean(gen,2,'omitnan')==9,:)=[];
time1=time; time1(mean(gen,2,'omitnan')==9,:)=[];
gen1=gen; gen1(mean(gen,2,'omitnan')==9,:)=[];
figure;
colormap(jet)
scatter(mean(dist1,2,'omitnan'),mean(time1,2,'omitnan'),10,mean(gen1,2,'omitnan'),'filled')
xlim([0 2000])
ylim([-25 110])
set(gca,'FontSize',14)
xlabel('Mean distance from index case (m)')
ylabel('Mean time since index case onset (months)')
cb=colorbar;
caxis([0 7])
title(cb,'Mean generation')
saveas(gcf,'MeanTimeVsMeanDistFromIndexCase.fig')
saveas(gcf,'MeanTimeVsMeanDistFromIndexCase.eps','epsc')
