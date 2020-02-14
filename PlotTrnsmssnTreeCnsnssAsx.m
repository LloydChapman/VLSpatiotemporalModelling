function PlotTrnsmssnTreeCnsnssAsx(rslts,nsmpls,infctn,infctr,src,infctnA,infctrA,srcA,iters,burnin1,varargin)

% N.B. Need to modify this to show consensus of which individuals had
% active KA at start of study

db='';
load(rslts)
clear M
% load(db)
load('data_final2.mat')

data=data(ismember(data.PARA,para),:);
% s=load(rslts);
% load(s.db)
% tmax=s.tmax;
% I1=s.I1;
% I=s.I;
% tI=s.tI;
% n=s.n;

if ~exist('z','var')
    z=burnin1+1:niters;
end

if nargin==10
    lon_rng=[min(data.HHNEWLNG) max(data.HHNEWLNG)];
    lat_rng=[min(data.HHNEWLAT) max(data.HHNEWLAT)];
else
    lon_rng=varargin{1};
    lat_rng=varargin{2};
end

lon=data.HHNEWLNG+(lon_rng(2)-lon_rng(1))/100*rand(n,1);
lat=data.HHNEWLAT+(lat_rng(2)-lat_rng(1))/100*rand(n,1);
% lon_rng=[min(lon) max(lon)];
% lat_rng=[min(lat) max(lat)];

scrnsz=get(0,'ScreenSize');

%% Most likely infectors and infection times
[Minfctr,Finfctr,Cinfctr]=mode(infctr,2);
Pinfctr=Finfctr/nsmpls;
[Minfctn,Finfctn,Cinfctn]=mode(infctn,2);
[Msrc,Fsrc,Csrc]=mode(src,2);

infctnA(infctnA==tmax+1)=NaN;
[MinfctrA,FinfctrA,CinfctrA]=mode(infctrA,2);
[MinfctnA,FinfctnA,CinfctnA]=mode(infctnA,2);
[MsrcA,FsrcA,CsrcA]=mode(src,2);

%% Transmission tree
% Set plotting colours
clrs=[[254 224 77]/255;[254 195 87]/255;[245 150 79]/255;0.8 0.255 0.145;[173 163 198]/255;[81 130 187]/255;[146 208 88]/255];
figure; 
set(gcf,'Units','Normalized','Position',[0 0 0.33 0.47])
set(0,'DefaultLegendAutoUpdate','off')
h=plot(lon,lat,'.','MarkerSize',8,'MarkerEdgeColor',clrs(1,:),'MarkerFaceColor',clrs(1,:)); hold on
% figure; h=plot(lon,lat,'.','MarkerSize',10); hold on
x_sclbar=lon_rng(1)+0.02*(lon_rng(2)-lon_rng(1));
y_sclbar=lat_rng(1)+0.02*(lat_rng(2)-lat_rng(1));
plot([x_sclbar,x_sclbar+0.001/(101/100)],[y_sclbar,y_sclbar],'k')
text(x_sclbar,y_sclbar+0.0002,'100m','FontSize',12)
set(gca,'FontSize',12,'XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[])
xlabel('Longitude'); ylabel('Latitude')
xlim(lon_rng); ylim(lat_rng)
i=(Minfctn<=maxIP);
% j=(Minfctn<=maxIP&Minfctr~=n+1);
ia=[];%(MinfctnA<=maxIP);
onset=tI(I);
InotP=setdiff(I,P);
trtmnt=tR(InotP);
trtmntIandP=tR(IandP);
% Add line here to overwrite treatment for VL case missing treatment time
trtmntAONR=mode(tRsAONR(:,z(iters)),2);
trtmntANONR=mode(tRsANONR(:,z(iters)),2);
onsetP=tP(I);
trtmntP=tRP(I);
k=(onset<=maxIP);
k1=(trtmnt<=maxIP);
k2=(trtmntIandP<=maxIP);
k3=(trtmntANONR<=maxIP);
k4=(trtmntAONR<=maxIP);
l=(onsetP<=maxIP);
l1=(trtmntP<=maxIP);
% hi=plot(lon(I1(i)),lat(I1(i)),'.','MarkerSize',14,'MarkerEdgeColor',[1 0.59 0],'MarkerFaceColor',[1 0.59 0]);
% % hia=plot(lon(ia),lat(ia),'y.','MarkerSize',14);
% hk=plot(lon(I(k)),lat(I(k)),'r.','MarkerSize',14);
% hk1=plot(lon(I(k1)),lat(I(k1)),'g.','MarkerSize',14);
% hl=plot(lon(I(l)),lat(I(l)),'m.','MarkerSize',14);
% hl1=plot(lon(I(l1)),lat(I(l1)),'g.','MarkerSize',14);

hi=plot(lon(I1(i)),lat(I1(i)),'.','MarkerSize',16,'MarkerEdgeColor',clrs(3,:),'MarkerFaceColor',clrs(3,:));
% hia=plot(lon(ia),lat(ia),'.','MarkerSize',16,'MarkerEdgeColor',clrs(2,:),'MarkerFaceColor',clrs(2,:));
hk=plot(lon(I(k)),lat(I(k)),'.','MarkerSize',16,'MarkerEdgeColor',clrs(4,:),'MarkerFaceColor',clrs(4,:));
% hk1=plot(lon(InotP(k1)),lat(InotP(k1)),'.','MarkerSize',16,'MarkerEdgeColor',clrs(7,:),'MarkerFaceColor',clrs(7,:));
plot(lon(InotP(k1)),lat(InotP(k1)),'.','MarkerSize',16,'MarkerEdgeColor',clrs(7,:),'MarkerFaceColor',clrs(7,:));
plot(lon(IandP(k2)),lat(IandP(k2)),'.','MarkerSize',16,'MarkerEdgeColor',clrs(5,:),'MarkerFaceColor',clrs(5,:));
% Should really distinguish initially potentially active VL cases who don't
% later develop PKDL from those that do, but there are only 2 cases of the latter.
plot(lon(ANONR(~k3)),lat(ANONR(~k3)),'.','MarkerSize',16,'MarkerEdgeColor',clrs(4,:),'MarkerFaceColor',clrs(4,:));
plot(lon(AONR(~k4)),lat(AONR(~k4)),'.','MarkerSize',16,'MarkerEdgeColor',clrs(4,:),'MarkerFaceColor',clrs(4,:));
plot(lon(ANONR(k3)),lat(ANONR(k3)),'.','MarkerSize',16,'MarkerEdgeColor',clrs(7,:),'MarkerFaceColor',clrs(7,:));
plot(lon(AONR(k4)),lat(AONR(k4)),'.','MarkerSize',16,'MarkerEdgeColor',clrs(7,:),'MarkerFaceColor',clrs(7,:));
hk1=plot(lon(setdiff(find(prevK),PI)),lat(setdiff(find(prevK),PI)),'.','MarkerSize',16,'MarkerEdgeColor',clrs(7,:),'MarkerFaceColor',clrs(7,:));
hk2=plot(lon(PI),lat(PI),'.','MarkerSize',16,'MarkerEdgeColor',clrs(5,:),'MarkerFaceColor',clrs(5,:));
hl=plot(lon(I(l)),lat(I(l)),'.','MarkerSize',16,'MarkerEdgeColor',clrs(6,:),'MarkerFaceColor',clrs(6,:));
hl1=plot(lon(I(l1)),lat(I(l1)),'.','MarkerSize',16,'MarkerEdgeColor',clrs(7,:),'MarkerFaceColor',clrs(7,:));
warning('off','MATLAB:legend:IgnoringExtraEntries');
% legend([h,hia,hi,hk,hk1,hk2,hl],'S','A','E','I','R','D','P')
legend([h,hi,hk,hk1,hk2,hl],'S/A','E','I','R','D','P')
% % ax = gca;
% % outerpos = ax.OuterPosition;
% % ti = ax.TightInset;
% % left = outerpos(1); %+ ti(1);
% % bottom = outerpos(2) + ti(2);
% % ax_width = outerpos(3)-0.2; % - ti(1) - ti(3);
% % ax_height = outerpos(4) - 1.5*ti(2);
% % ax.Position = [left bottom ax_width ax_height];
% % set(gcf,'Units','normalized')
% % set(gcf,'Position',[0.305 0.419 0.32 0.47])
vid = VideoWriter('ExmplTrnsmssnTreeCnsnss.avi');
% vid = VideoWriter('ExmplTrnsmssnTreeCnsnss.mp4','MPEG-4'); % only works on Mac/Windows
vid.FrameRate=2;
open(vid)
for t=maxIP+1:tmax
%     j=(infctntime==i);
    i=(Minfctn==t);
    ia=(MinfctnA==t&MinfctrA~=n+1);
    j=find(Minfctn==t&Minfctr~=n+1);
    k=(onset==t);
    k1=(trtmnt==t);
    k2=(trtmntIandP==t);
    l=(onsetP==t); 
    l1=(trtmntP==t);
%     if t==1 || any(i)
%     hi=plot(lon(I1(i)),lat(I1(i)),'.','MarkerSize',14,'MarkerEdgeColor',[1 0.59 0],'MarkerFaceColor',[1 0.59 0]);
%     end
% %     if t==1 || any(ia)
% %     hia=plot(lon(ia),lat(ia),'y.','MarkerSize',14);
% %     end
%     if t==1 || any(k)
%     hk=plot(lon(I(k)),lat(I(k)),'r.','MarkerSize',14);
%     end
%     if t==1 || any(k1)
%     hk1=plot(lon(I(k1)),lat(I(k1)),'g.','MarkerSize',14);
%     end
%     if t==1 || any(l)
%     hl=plot(lon(I(l)),lat(I(l)),'m.','MarkerSize',32);
%     end
%     if t==1 || any(l1)
%     hl1=plot(lon(I(l1)),lat(I(l1)),'g.','MarkerSize',14);
%     end

%     if t==1 || any(ia)
%     hia=plot(lon(ia),lat(ia),'.','MarkerSize',14,'MarkerEdgeColor',clrs(2,:),'MarkerFaceColor',clrs(2,:));
%     end
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
        quiver_thick(lon(Minfctr(m)),lat(Minfctr(m)),0.95*(lon(I1(m))-lon(Minfctr(m))),0.95*(lat(I1(m))-lat(Minfctr(m))),'arrow_thickness',1e-4/3,'arrow_thickness_option','absolute','arrowhead_length',2e-4,'arrowhead_length_option','absolute','axis_equal',1,'plot_color',clrs(Msrc(m),:),'EdgeColor',clrs(Msrc(m),:),'EdgeAlpha',0,'FaceAlpha',Pinfctr(m)^(2/3),'LineWidth',1); 
    end
    xlim(lon_rng); ylim(lat_rng)
%     if any(i) && any(k) && any(k1) && any(l)
%         legend([h,hia,hi,hk,hk1,hk2,hl],'S','A','E','I','R','D','P')
        legend([h,hi,hk,hk1,hk2,hl],'S/A','E','I','R','D','P')
        title(['t = ' num2str(t)])
%     end
    pause(0.5)
    if ismember(t,[24,48,96])
%         saveas(gcf,['TrnsmssnTreeMonth' num2str(t) 'Para1.png'])
        saveas(gcf,['ExmplTrnsmssnTreeMonth' num2str(t)])
        saveas(gcf,['ExmplTrnsmssnTreeMonth' num2str(t) '.eps'],'epsc')
    end
    M(t)=getframe(gcf);
    writeVideo(vid,M(t));
end
hold off
close(vid)

% %% Infection sources
% % Frequencies of different infection sources
% edges=[2,3,4,6,9,Inf];
% srcfreq=NaN(numel(I1),numel(edges)-1);
% for i=1:numel(I1)
%     srcfreq(i,:)=histcounts(src(i,:),edges);
% end
% srcfreq=srcfreq/nsmpls;
% % Stacked bar chart showing proportion of samples for which most likely
% % infector was from given source (bckgrnd, VL, PKDL)
% [~,ord]=sort(tI(I1));
% srcfreqsrt=srcfreq(ord,[5,1:4]);
% figure; 
% set(gcf, 'Position', [0 70 scrnsz(3) scrnsz(4) - 150]);
% h1=bar(srcfreqsrt,1,'stacked');
% h1(1).FaceColor=[0.96 0.9 0.8]; %[0.965 0.698 0.4196];
% h1(1).FaceAlpha=0.9;
% h1(2).FaceColor=[254 195 87]/255;
% h1(2).FaceAlpha=0.8;
% h1(3).FaceColor=[245 150 79]/255;
% h1(3).FaceAlpha=0.8;
% h1(4).FaceColor=[0.8 0.255 0.145];
% h1(4).FaceAlpha=0.8;
% h1(5).FaceColor=[81 130 187]/255;
% set(gca,'FontSize',20)
% xlabel('Case number (by onset)'); ylabel('Probability of infection source')
% legend('Bckgrnd','Asx','Presx','VL','PKDL','Location','northwest')
% xlim([0.5 numel(I1)+0.5])
% ylim([0 1])
% % saveas(gcf,'ProbInfctnSourceAsx.png','png')
% 
% %% Reproduction numbers
% % VL
% infctrfrqs=zeros(nI,nsmpls);
% for i=1:size(infctr,2)
%     tmp=histcounts(infctr(src(:,i)==4,i),1:n+1);
%     infctrfrqs(:,i)=tmp(I);
% end
% infctrAfrqs=zeros(nI,nsmpls);
% for i=1:size(infctrA,2)
%     tmp=histcounts(infctrA(srcA(:,i)==4,i),1:n+1);
%     infctrAfrqs(:,i)=tmp(I);
% end
% Ri=infctrfrqs+infctrAfrqs;
% meanRi=mean(Ri,2);
% HPDIRi=zeros(nI,2);
% for i=1:nI
% [~,HPDIRi(i,:)]=CalcModeAndHPDI(Ri(i,:),[],'auto');
% end
% infctrfreq=histcounts(infctr(src==4),1:n+1);
% infctrfreq=infctrfreq(I)/nsmpls; % exclude potentially active KA cases at start of study for now as they should only be counted if they had active KA during the study period 
% % infctrfreq=infctrfreq([I;A])/nsmpls;
% % Reproduction numbers of cases
% figure; bar(infctrfreq)
% xlabel('VL case number'); ylabel('Reproduction number')
% % Distribution of reproduction numbers of cases
% figure; histogram(infctrfreq,0:0.5:max(infctrfreq)+1,'Normalization','pdf')
% xlim([0 Inf])
% set(gca,'FontSize',14)
% xlabel('Case reproduction number'); ylabel('Density')
% saveas(gcf,'VLReffDistn.png','png')
% % mean(infctrfreq)
% % median(infctrfreq)
% % prctile(infctrfreq,[2.5 97.5])
% % Sort reproduction numbers by onset month
% [~,ord1]=sort(tI(I));
% % [~,ord1]=sort(tI([I;A]));
% infctrfreqsrt=infctrfreq(ord1);
% % Case reproduction numbers over time
% figure; bar(infctrfreqsrt);
% xlabel('VL case number (by onset)'); ylabel('Case reproduction number')
% 
% % PKDL
% PKDLinfctrfreq=histcounts(infctr(src==6),1:n+1);
% PKDLinfctrfreq=PKDLinfctrfreq(P)/nsmpls;
% % Reproduction numbers of cases
% figure; bar(PKDLinfctrfreq)
% xlabel('PKDL case number'); ylabel('Reproduction number')
% % Distribution of reproduction numbers of cases
% figure; histogram(PKDLinfctrfreq,'Normalization','pdf')
% xlim([0 Inf])
% set(gca,'FontSize',14)
% xlabel('Case reproduction number'); ylabel('Density')
% saveas(gcf,'PKDLReffDistn.png','png')
% % mean(infctrfreq)
% % median(infctrfreq)
% % prctile(infctrfreq,[2.5 97.5])
% % Sort reproduction numbers by onset month
% [~,ord2]=sort(tP(P));
% PKDLinfctrfreqsrt=PKDLinfctrfreq(ord2);
% % Case reproduction numbers over time
% figure; bar(PKDLinfctrfreqsrt);
% xlabel('PKDL case number (by onset)'); ylabel('Reproduction number')
% 
% %% Infection distances
% % Over all samples and time
% figure; histogram(dist,50,'Normalization','pdf')
% xlim([0 300])
% set(gca,'FontSize',14)
% xlabel('Infection distance (m)'); ylabel('Density')
% saveas(gcf,'InfctnDstnceDistn.png','png')
% % In each year
% figure; xlabel('Distance (m)'); ylabel('Density')
% hold on; 
% for i=1:ceil(max(max(infctn))/12)
% histogram(dist(infctn>(i-1)*12&infctn<=i*12),'Normalization','pdf'); pause(1)
% end
% hold off
% 
% %% Times between onsets and secondary infections
% % VL
% figure; histogram(times(src==3),'Normalization','pdf')
% xlabel('Time between VL onset and secondary infection (months)'); ylabel('Density')
% 
% % PKDL
% figure; histogram(times(src==5),'Normalization','pdf')
% xlabel('Time between PKDL onset and secondary infection (months)'); ylabel('Density')
% 
% %% Relationship between infection distances and times to secondary infections
% % VL
% figure; histogram2(dist(src==3),times(src==3),'Normalization','pdf','DisplayStyle','tile'); hcb1=colorbar; title(hcb1,'pmf')
% xlabel('Infection distance (m)'); ylabel('Time between VL onset and secondary infection (months)')
% 
% % PKDL
% figure; histogram2(dist(src==5),times(src==5),'Normalization','pdf','DisplayStyle','tile'); hcb2=colorbar; title(hcb2,'pmf')
% xlabel('Infection distance (m)'); ylabel('Time between PKDL onset and secondary infection (months)')
% 
% %% Relationship between VL/PKDL onset-to-treatment and infection distances and times
% OT=trtmntinfctr-onsetinfctr;
% 
% % VL
% figure; histogram2(OT(src==3),times(src==3),'Normalization','pdf','DisplayStyle','tile'); hcb1=colorbar; title(hcb1,'pmf')
% xlabel('VL onset-to-treatment time (months)'); ylabel('Time between VL onset and secondary infection (months)')
% 
% % PKDL
% figure; histogram2(OT(src==5),times(src==5),'Normalization','pdf','DisplayStyle','tile'); hcb2=colorbar; title(hcb2,'pmf')
% xlabel('PKDL onset-to-resolution time (months)'); ylabel('Time between PKDL onset and secondary infection (months)')
% 
% %% VL/PKDL onset-to-treatment time vs reproduction number
% OTcase=tR(I)-tI(I);
% x=OTcase(~isnan(OTcase));
% y=infctrfreq(~isnan(OTcase))';
% b1=x\y;
% rsq=1-sum((y-b1*x).^2)/sum((y-mean(y)).^2);
% figure; plot(tR(I)-tI(I),infctrfreq,'.','MarkerSize',12); hold on
% plot(x,b1*x,'LineWidth',2)
% set(gca,'FontSize',14)
% xlabel('VL infector onset-to-treatment time (months)'); ylabel('Case reproduction number')
% % saveas(gcf,'VLinfctrOTvsReff.png','png')
% 
% figure; plot(tRP(P)-tP(P),PKDLinfctrfreq,'.',tmax-tP(~isnan(tP)&isnan(tRP)&~PothrObs),PKDLinfctrfreq(isnan(tRP(P))),'.','MarkerSize',12)
% set(gca,'FontSize',14)
% xlabel('Infector onset-to-treatment time (months)'); ylabel('Case reproduction number')
% legend('Resolved','Unresolved')
% % saveas(gcf,'PKDLinfctrOTvsReff.png','png')
