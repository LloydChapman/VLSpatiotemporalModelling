function RunPlotDevDistn(rslts,burnin,IPD,str,h0s,h40s,delta0s)
figure;
nMdls=numel(rslts);

lgd=cell(1,nMdls);
for i=1:nMdls
    hf(i)=PlotDevDistn(rslts{i},burnin); hold on
    if delta0s(i)==0
        str1='$$\delta=0, ';
    else
        str1='$$\delta>0, ';
    end        
    lgd{i}=[str1 'h_0=' num2str(h0s(i)) ', h_4=' num2str(h40s(i)) '$$'];
end
% legend(hL,cellfun(@(x)x(11:end),rslts,'UniformOutput',false),'Interpreter','none')
legend(hf,lgd,'Interpreter','LaTeX')
saveas(gcf,['PstrDevDistns' IPD str])
saveas(gcf,['PstrDevDistns' IPD str '.eps'],'epsc')
saveaspdf(gcf,['PstrDevDistns' IPD str])