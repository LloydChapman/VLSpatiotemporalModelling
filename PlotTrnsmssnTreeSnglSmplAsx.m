function [psmpl,p1smpl]=PlotTrnsmssnTreeSnglSmplAsx(rslts,infctn,infctr,src,infctnA,infctrA,srcA,iters,burnin1,varargin)

load(rslts)
clear M
load('data_final2.mat')

data=data(ismember(data.PARA,para),:);

if ~exist('z','var')
    z=burnin1+1:niters;
end

% Calculate distances of parameter sets from modes
nbins=50;
mode_p=NaN(1,np);
for i=1:np
    mode_p(i)=CalcModeAndHPDI(p(z(iters),i),nbins);
end
mode_p1=CalcModeAndHPDI(p1(z(iters)),nbins);
x=sum(abs(bsxfun(@minus,[p(z(iters),:),p1(z(iters))],[mode_p,mode_p1])./[mode_p,mode_p1]),2);

if nargin==9
    lon_rng=[min(data.HHNEWLNG) max(data.HHNEWLNG)];
    lat_rng=[min(data.HHNEWLAT) max(data.HHNEWLAT)];
    % Find sample for which parameter values are closest to modes
    smpl=find(x==min(x),1);
elseif nargin==10
    lon_rng=[min(data.HHNEWLNG) max(data.HHNEWLNG)];
    lat_rng=[min(data.HHNEWLAT) max(data.HHNEWLAT)];
    smpl=varargin{1};
elseif nargin==11
    lon_rng=varargin{1};
    lat_rng=varargin{2};
    smpl=find(x==min(x),1);
elseif nargin==12
    lon_rng=varargin{1};
    lat_rng=varargin{2};
    smpl=varargin{3};
end
psmpl=p(z(iters(smpl)),:);
p1smpl=p1(z(iters(smpl)));

lon=data.HHNEWLNG+(lon_rng(2)-lon_rng(1))/100*rand(n,1);
lat=data.HHNEWLAT+(lat_rng(2)-lat_rng(1))/100*rand(n,1);

scrnsz=get(0,'ScreenSize');

%% Extract infectors and infection times for chosen sample
infctr=infctr(:,smpl);
infctn=infctn(:,smpl);
src=src(:,smpl);

infctrA=infctrA(:,smpl);
infctnA=infctnA(:,smpl);
srcA=srcA(:,smpl);

%% Transmission tree
% Set plotting colours
clrs=[[254 224 77]/255;[254 195 87]/255;[245 150 79]/255;0.8 0.255 0.145;[173 163 198]/255;[81 130 187]/255;[146 208 88]/255];
figure; %set(gcf,'Units','Normalized','Position',[0 0 0.5 0.8]) % uncomment to plot larger
set(gcf,'Units','Normalized','Position',[0 0 0.37125 0.52875])
set(0,'DefaultLegendAutoUpdate','off')
h=plot(lon,lat,'.','MarkerSize',8,'MarkerEdgeColor',clrs(1,:),'MarkerFaceColor',clrs(1,:)); hold on
x_sclbar=lon_rng(1)+0.02*(lon_rng(2)-lon_rng(1));
y_sclbar=lat_rng(1)+0.02*(lat_rng(2)-lat_rng(1));
plot([x_sclbar,x_sclbar+0.001/(101/100)],[y_sclbar,y_sclbar],'k')
text(x_sclbar,y_sclbar+0.0002,'100m','FontSize',12)
set(gca,'FontSize',12,'XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[])
xlabel('Longitude'); ylabel('Latitude')
xlim(lon_rng); ylim(lat_rng)
i=(infctn<=maxIP);
ia=(infctnA<=maxIP);
onset=tI(I);
InotP=setdiff(I,P);
trtmntA=tRAs(:,z(iters(smpl)));
trtmnt=tR(InotP);
trtmntIandP=tR(IandP);
% Add line here to overwrite treatment for VL case missing treatment time
trtmntAONR=tRsAONR(:,z(iters(smpl)));
trtmntANONR=tRsANONR(:,z(iters(smpl)));
onsetP=tP(I);
trtmntP=tRP(I);
ra=(trtmntA<=maxIP);
k=(onset<=maxIP);
k1=(trtmnt<=maxIP);
k2=(trtmntIandP<=maxIP);
k3=(trtmntANONR<=maxIP);
k4=(trtmntAONR<=maxIP);
l=(onsetP<=maxIP);
l1=(trtmntP<=maxIP);

hia=plot(lon(ia),lat(ia),'.','MarkerSize',8,'MarkerEdgeColor',clrs(2,:),'MarkerFaceColor',clrs(2,:));
hra=plot(lon(ra),lat(ra),'.','MarkerSize',8,'MarkerEdgeColor',clrs(7,:),'MarkerFaceColor',clrs(7,:));
hi=plot(lon(I1(i)),lat(I1(i)),'.','MarkerSize',16,'MarkerEdgeColor',clrs(3,:),'MarkerFaceColor',clrs(3,:));
hk=plot(lon(I(k)),lat(I(k)),'.','MarkerSize',16,'MarkerEdgeColor',clrs(4,:),'MarkerFaceColor',clrs(4,:));
plot(lon(InotP(k1)),lat(InotP(k1)),'.','MarkerSize',16,'MarkerEdgeColor',clrs(7,:),'MarkerFaceColor',clrs(7,:));
plot(lon(IandP(k2)),lat(IandP(k2)),'.','MarkerSize',16,'MarkerEdgeColor',clrs(5,:),'MarkerFaceColor',clrs(5,:));
plot(lon(ANONR(~k3)),lat(ANONR(~k3)),'.','MarkerSize',16,'MarkerEdgeColor',clrs(4,:),'MarkerFaceColor',clrs(4,:)); % should really distinguish initially potentially active VL cases who don't later develop PKDL from those that do, but there are only 2 cases of the latter
plot(lon(AONR(~k4)),lat(AONR(~k4)),'.','MarkerSize',16,'MarkerEdgeColor',clrs(4,:),'MarkerFaceColor',clrs(4,:));
plot(lon(ANONR(k3)),lat(ANONR(k3)),'.','MarkerSize',16,'MarkerEdgeColor',clrs(7,:),'MarkerFaceColor',clrs(7,:));
plot(lon(AONR(k4)),lat(AONR(k4)),'.','MarkerSize',16,'MarkerEdgeColor',clrs(7,:),'MarkerFaceColor',clrs(7,:));
hk1=plot(lon(setdiff(find(prevK),PI)),lat(setdiff(find(prevK),PI)),'.','MarkerSize',16,'MarkerEdgeColor',clrs(7,:),'MarkerFaceColor',clrs(7,:));
hk2=plot(lon(PI),lat(PI),'.','MarkerSize',16,'MarkerEdgeColor',clrs(5,:),'MarkerFaceColor',clrs(5,:));
hl=plot(lon(I(l)),lat(I(l)),'.','MarkerSize',16,'MarkerEdgeColor',clrs(6,:),'MarkerFaceColor',clrs(6,:));
hl1=plot(lon(I(l1)),lat(I(l1)),'.','MarkerSize',16,'MarkerEdgeColor',clrs(7,:),'MarkerFaceColor',clrs(7,:));
warning('off','MATLAB:legend:IgnoringExtraEntries');
legend([h,hia,hi,hk,hk1,hk2,hl],'S','A','E','I','R','D','P')
% legend([h,hi,hk,hk1,hk2,hl],'S/A','E','I','R','D','P')
vid = VideoWriter('ExmplTrnsmssnTreeSnglSmpl.avi');
% vid = VideoWriter('ExmplTrnsmssnTreeSnglSmpl.mp4','MPEG-4'); % only works on Mac/Windows
vid.FrameRate=2;
open(vid)
for t=maxIP+1:tmax
    i=(infctn==t);
    ia=(infctnA==t);
    ra=(trtmntA==t);
    j=find(infctn==t&infctr~=n+1);
    k=(onset==t);
    k1=(trtmnt==t);
    k2=(trtmntIandP==t);
    l=(onsetP==t); 
    l1=(trtmntP==t);

    if t==1 || any(ia)
    hia=plot(lon(ia),lat(ia),'+','MarkerSize',10,'MarkerEdgeColor',clrs(2,:),'MarkerFaceColor',clrs(2,:));
    end
    if t==1 || any(ra)
    hra=plot(lon(ra),lat(ra),'+','MarkerSize',10,'MarkerEdgeColor',clrs(7,:),'MarkerFaceColor',clrs(7,:));
    end
    if t==1 || any(i)
    hi=plot(lon(I1(i)),lat(I1(i)),'.','MarkerSize',16,'MarkerEdgeColor',clrs(3,:),'MarkerFaceColor',clrs(3,:));
    end
    if t==1 || any(k)
    hk=plot(lon(I(k)),lat(I(k)),'.','MarkerSize',16,'MarkerEdgeColor',clrs(4,:),'MarkerFaceColor',clrs(4,:));
    end
    if t==1 || any(k1)
    hk1=plot(lon(InotP(k1)),lat(InotP(k1)),'.','MarkerSize',16,'MarkerEdgeColor',clrs(7,:),'MarkerFaceColor',clrs(7,:));
    end
    if t==1 || any(k2)
    hk2=plot(lon(IandP(k2)),lat(IandP(k2)),'.','MarkerSize',16,'MarkerEdgeColor',clrs(5,:),'MarkerFaceColor',clrs(5,:));
    end
    if t==1 || any(l)
    hl=plot(lon(I(l)),lat(I(l)),'.','MarkerSize',16,'MarkerEdgeColor',clrs(6,:),'MarkerFaceColor',clrs(6,:));
    end
    if t==1 || any(l1)
    hl1=plot(lon(I(l1)),lat(I(l1)),'.','MarkerSize',16,'MarkerEdgeColor',clrs(7,:),'MarkerFaceColor',clrs(7,:));
    end
    for ii=1:numel(j)
        m=j(ii);
        % Scale x and y components of arrow between infector and infectee so that arrow head is visible
%         quiver(lon(Minfctr(m)),lat(Minfctr(m)),0.95*(lon(I1(m))-lon(Minfctr(m))),0.95*(lat(I1(m))-lat(Minfctr(m))),0,'Color',clrs(Msrc(m),:),'LineWidth',1,'MaxHeadSize',0.4)
        % Options for plotting at default MATLAB size
        quiver_thick(lon(infctr(m)),lat(infctr(m)),0.95*(lon(I1(m))-lon(infctr(m))),0.95*(lat(I1(m))-lat(infctr(m))),'arrow_thickness',1e-4/3,'arrow_thickness_option','absolute','arrowhead_length',2e-4,'arrowhead_length_option','absolute','axis_equal',1,'plot_color',clrs(src(m),:),'EdgeColor',clrs(src(m),:),'EdgeAlpha',0,'LineWidth',1); 
        % Options for plotting at larger size
%         quiver_thick(lon(infctr(m)),lat(infctr(m)),0.99*(lon(I1(m))-lon(infctr(m))),0.99*(lat(I1(m))-lat(infctr(m))),'arrow_thickness',1e-4/6,'arrow_thickness_option','absolute','arrowhead_length',1e-4,'arrowhead_length_option','absolute','axis_equal',1,'plot_color',clrs(src(m),:),'EdgeColor',clrs(src(m),:),'EdgeAlpha',0,'LineWidth',1);
    end
    xlim(lon_rng); ylim(lat_rng)
%     if any(i) && any(k) && any(k1) && any(l)
        legend([h,hia,hi,hk,hk1,hk2,hl],'S','A','E','I','R','D','P')
%         legend([h,hi,hk,hk1,hk2,hl],'S/A','E','I','R','D','P')
        title(['t = ' num2str(t)])
%     end
    pause(0.5)
    if ismember(t,[24,48,96])
        saveas(gcf,['ExmplTrnsmssnTreeSnglSmplMonth' num2str(t)])
        saveas(gcf,['ExmplTrnsmssnTreeSnglSmplMonth' num2str(t) '.eps'],'epsc')
    end
    M(t)=getframe(gcf);
    writeVideo(vid,M(t));
end
hold off
close(vid)
