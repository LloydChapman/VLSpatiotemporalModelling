function PlotMeanIndexCaseToInfcteeDistsAndTimes(rslts,gen,dist,time,chain)

load(rslts,'n')
nsmpls=size(gen,2);
idx=NaN(size(chain,1),nsmpls);
for j=1:nsmpls
    for i=1:size(chain,1)
        if ~isempty(chain{i,j})
            idx(i,j)=chain{i,j}(end);
        end
    end
end

%% Calculate mean index-case-to-infectee distances and index-case-onset-to-infectee-infection times for each index case
ngen=max(max(gen));
unqinfctr=setdiff(unique(idx(~isnan(idx))),n+1);
ninfctr=numel(unqinfctr);
meandist=NaN(ninfctr,ngen,nsmpls);
meantime=NaN(ninfctr,ngen,nsmpls);
for k=1:nsmpls
    for i=1:ninfctr
        disti=dist(idx(:,k)==unqinfctr(i),k);
        timei=time(idx(:,k)==unqinfctr(i),k);
        geni=gen(idx(:,k)==unqinfctr(i),k);
        for j=1:max(geni)
            meandist(i,j,k)=mean(disti(geni==j));
            meantime(i,j,k)=mean(timei(geni==j));
        end
    end
end

% Calculate mean index-case-to-infectee distances and times in each 
% generation over all sampled transmission trees
meanmeandist=mean(meandist,3,'omitnan');
meanmeantime=mean(meantime,3,'omitnan');

%% Plot mean index-case-onset-to-infectee-infection times vs mean index-case-to-infectee distances
figure; hold on
set(gca,'FontSize',14)
cmap=jet(ngen);
lgd=cell(1,ngen);
for j=1:ngen
    plot(meanmeandist(:,j),meanmeantime(:,j),'.','Color',cmap(ngen+1-j,:),'MarkerSize',12)
    lgd{j}=num2str(j);
end
hold off
lh=legend(lgd);
title(lh,'Generation','FontSize',10)
xlabel('Mean index-case-to-infectee distance (m)')
ylabel('Mean index-case-onset-to-infectee-infection time (months)')
xlim([0 1000])
ylim([0 max(max(meanmeantime))+1])
saveas(gcf,'MeanIndexCaseToInfcteeTimeVsDist')
saveas(gcf,'MeanIndexCaseToInfcteeTimeVsDist.eps','epsc')