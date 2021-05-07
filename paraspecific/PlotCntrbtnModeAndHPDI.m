function PlotCntrbtnModeAndHPDI(x,nbins,nsrcs,t,clrs,xlbl,ylbl,mode,lgd,varargin)
nrow=size(x,2);
M=zeros(nrow,nsrcs);
HPDI=zeros(nrow,2,nsrcs);
for i=1:nsrcs
    for j=1:nrow
        if mode
            [M(j,i),HPDI(j,:,i)]=CalcModeAndHPDI(x(:,j,i),nbins);
        else
            M(j,i)=median(x(:,j,i));
            HPDI(j,:,i)=[quantile(x(:,j,i),0.025) quantile(x(:,j,i),0.975)];
        end
    end
end
t2=[t,fliplr(t)];
figure; hold on
for i=1:nsrcs
hf1(i)=plot(t,M(:,i),'LineWidth',1,'Color',clrs(i,:));
HPDR=[HPDI(:,1,i)',fliplr(HPDI(:,2,i)')];
hf2=fill(t2,HPDR,clrs(i,:),'LineStyle','none');
hf2.FaceAlpha=0.5;
end
set(gca,'FontSize',14)
xlabel(xlbl); ylabel(ylbl)
xlim([t(1) t(end)])
if nargin==10
    lcn=varargin{1};
else
    lcn='northeast';
end
legend(hf1,lgd,'Location',lcn)