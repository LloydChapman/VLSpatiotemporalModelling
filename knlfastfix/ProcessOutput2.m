function [mode_p,HPDI,mode_p1,HPDI1,pcorr,pIPtAcorr,mode_sptl,HPDI_sptl,mode_bckgrnd,HPDI_bckgrnd,mode_d_half,HPDI_d_half,mode_d_half_out,HPDI_d_half_out,mode_WHHRI,HPDI_WHHRI,mean_IP,mode_IP,HPDI_IP]=ProcessOutput2(str,burnin1,thin,plotMltpl,doPlots,savePlots)
close all

% Load MCMC output
load(str)
rslts=rslts(6:end); % remove "MCMC_" from rslts string

% Discard burnins
if nargin==1
    z=burnin+1:niters; % if no burnin supplied, use saved default (round(niters/10))
    thin=1;
    plotMltpl=true;
    doPlots=false;
    savePlots=false;
else
    z=burnin1+1:niters; % if new burnin supplied, use it
end

% Load data
load('data_final2.mat')
 
% Select data for para
data=data(ismember(data.PARA,para),:);
 
% Rename longitude and latitude variables
data.Properties.VariableNames{'HHNEWLNG'}='longitude';
data.Properties.VariableNames{'HHNEWLAT'}='latitude';

%% TRACE PLOTS AND POSTERIOR DISTRIBUTIONS OF LOG-LIKELIHOOD AND PARAMETERS
nbins=50;
scrnsz=get(0,'ScreenSize');
zthin=z(1:thin:end);
nz=numel(zthin);
% Plot output with fixed parameters excluded
figure;
if plotMltpl % multiple realisations of missing KA onset, treatment, relapse and relapse treatment times
    [mode_pu,HPDIu,mode_p1u,HPDI1u]=PlotOutput2(zthin,LL,p,nu,pname,priorp,p1,a,b,n,tmax,I,RpreD,DpreR,OR,NONR,RNO,ONR,A,ANONR,AONR,RL,RLO,RLNO,tI,tR,tD,tRL,tRLR,tIsNONR,tIsRNO,tRsNONR,tRsONR,tIsANONR,tRsANONR,tRsAONR,tRLsRLNO,tRLRsRLO,tRLRsRLNO,nbins,scrnsz);
else % 1 realisation of missing times
    [mode_p,HPDI,mode_p1,HPDI1]=PlotOutput(zthin,LL,p,nu,pname,priorp,p1,a,b,n,tmax,I,RpreD,DpreR,RL,tI,tR,tD,tRL,tRLR,nbins,scrnsz);
end
saveas2(gcf,['PSTR_DISTNS_' rslts],savePlots)
saveas2(gcf,['PSTR_DISTNS_' rslts '.eps'],savePlots,'epsc')
figure;
PlotTrace(zthin,p,nu,pname,p1,mode_pu,HPDIu,mode_p1u,HPDI1u,scrnsz)

% Plot output with fixed parameters included
figure;
if plotMltpl % multiple realisations of missing times
    [mode_p,HPDI,mode_p1,HPDI1]=PlotOutput2(zthin,LL,p,np,pname,priorp,p1,a,b,n,tmax,I,RpreD,DpreR,OR,NONR,RNO,ONR,A,ANONR,AONR,RL,RLO,RLNO,tI,tR,tD,tRL,tRLR,tIsNONR,tIsRNO,tRsNONR,tRsONR,tIsANONR,tRsANONR,tRsAONR,tRLsRLNO,tRLRsRLO,tRLRsRLNO,nbins,scrnsz);
else % 1 realisation of missing times
    [mode_p,HPDI,mode_p1,HPDI1]=PlotOutput(zthin,LL,p,np,pname,priorp,p1,a,b,n,tmax,I,RpreD,DpreR,RL,tI,tR,tD,tRL,tRLR,nbins,scrnsz);
end
figure;
PlotTrace(zthin,p,np,pname,p1,mode_p,HPDI,mode_p1,HPDI1,scrnsz)

%% AUTOCORRELATION AND PARAMETER CORRELATION
% Autocorrelation plots for each parameter
figure; set(gcf, 'Position', [0 70 round(scrnsz(3)/2) round(scrnsz(3)/2.5)]);
for i=1:nu
    j=u(i);
    subplot(ceil((nu+1)/2),2,i)
    acf(p(zthin,j),min(2000,nz-1));
    title(['\' pname{j}])
end
subplot(ceil((nu+1)/2),2,nu+1)
% figure;
acf(p1(zthin),min(2000,nz-1));
title('$$p$$','Interpreter','latex')
saveas2(gcf,['ACFs_' rslts],savePlots)
saveas2(gcf,['ACFs_' rslts '.eps'],savePlots,'epsc')

% Matrix of parameter pair plots
figure; set(gcf, 'Position', [0 70 round(scrnsz(3)/2) round(scrnsz(3)/2.5)]);
[~,ax,~,~,~]=plotmatrix2([p(zthin,u),p1(zthin)]);
for i=1:nu
    j=u(i);
    xlabel(ax(end,i),['\' pname{j}],'Fontsize',13)
    ylabel(ax(i,1),['\' pname{j} '    '],'Fontsize',13,'rot',0)
end
xlabel(ax(end,nu+1),'$$p$$','Fontsize',13,'Interpreter','latex')
ylabel(ax(nu+1,1),'$$p$$    ','Fontsize',13,'Interpreter','latex','rot',0)
saveas2(gcf,['ParamCrrltn_' rslts],savePlots)
saveas2(gcf,['ParamCrrltn_' rslts '.eps'],savePlots,'epsc')

% Parameter correlation coefficients
pcorr=corrcoef([p(zthin,u) p1(zthin)]);

%% CORRELATION BETWEEN TRANSMISSION PARAMETERS AND ASYMPTOMATIC INFECTION TIMES
% Correlation between beta and mean asymptomatic infection time
mean_tA=NaN(nz,1);
for i=1:nz
    mean_tA(i)=mean(tAs(tAs(:,zthin(i))>0 & tAs(:,zthin(i))<tmax+1,zthin(i)));
end
figure; [~,ax,~,~,~]=plotmatrix([p(zthin,1),mean_tA]);
xlabel(ax(end,1),'\beta','Fontsize',13)
ylabel(ax(1,1),'\beta    ','Fontsize',13,'rot',0)
xlabel(ax(end,2),'$$\bar{A}$$','Fontsize',13,'Interpreter','latex')
ylabel(ax(2,1),'$$\bar{A}$$    ','Fontsize',13,'Interpreter','latex','rot',0)
saveas2(gcf,['betaAsxInfctnTimeCrrltn' rslts],savePlots)
saveas2(gcf,['betaAsxInfctnTimeCrrltn' rslts '.eps'],savePlots,'epsc')

% Correlation between transmission parameters, incubation period
% distribution parameter, incubation periods and asymptomatic infection
% times
pIPtA=[p(zthin,u),p1(zthin),mean(IPs(:,zthin),1)',mean_tA];
figure; set(gcf, 'Position', [0 70 round(scrnsz(3)/2) round(scrnsz(3)/2.5)]);
[~,ax,~,~,~]=plotmatrix(pIPtA);
for i=1:nu
    j=u(i);
    xlabel(ax(end,i),['\' pname{j}],'Fontsize',13)
    ylabel(ax(i,1),['\' pname{j} '    '],'Fontsize',13,'rot',0)
end
xlabel(ax(end,nu+1),'$$p$$','Fontsize',13,'Interpreter','latex')
ylabel(ax(nu+1,1),'$$p$$    ','Fontsize',13,'Interpreter','latex','rot',0)
xlabel(ax(end,nu+2),'$$\bar{IP}$$','Fontsize',13,'Interpreter','latex')
ylabel(ax(nu+2,1),'$$\bar{IP}$$    ','Fontsize',13,'Interpreter','latex','rot',0)
xlabel(ax(end,nu+3),'$$\bar{A}$$','Fontsize',13,'Interpreter','latex')
ylabel(ax(nu+3,1),'$$\bar{A}$$    ','Fontsize',13,'Interpreter','latex','rot',0)
pIPtAcorr=corrcoef(pIPtA);
saveas2(gcf,['pIPtACrrltn' rslts],savePlots)
saveas2(gcf,['pIPtACrrltn' rslts '.eps'],savePlots,'epsc')

%% ESTIMATED TRANSMISSION KERNEL
% Plot estimated transmission kernel
dHH=CalcHHDists(data);
if strcmp(typ,'Cauchy')
    dHHsqrd=dHH.^2;
else
    dHHsqrd=[];
end
[mode_Ke,mode_K0,mode_rate,HPDI_rate]=PlotKnl2(zthin,p,K0,mode_p,HPDI,dHH,dHHsqrd,typ,n,nHH,f);

% Save plot
if ~all(isnan(mode_Ke))
    saveas2(gcf,['SPTL_KRNL_' rslts],savePlots)
    saveas2(gcf,['SPTL_KRNL_' rslts '.eps'],savePlots,'epsc')
    % Remove within-HH part of spatial kernel
    xl=xlim;
    xlim([1.1 xl(2)])
    legend('hide')
    saveas2(gcf,['SPTL_KRNL_' rslts '_OUTSIDE_HH'],savePlots)
    saveas2(gcf,['SPTL_KRNL_' rslts '_OUTSIDE_HH'],savePlots,'epsc')
end

% Calculate transmission rate from VL cases and background transmission rate
mode_sptl=mode_rate(1)*1e4; % cases/10,000 people/mnth
HPDI_sptl=HPDI_rate*1e4; % cases/10,000 people/mnth
mode_bckgrnd=mode_p(3)*1e4; % cases/10,000 people/mnth
HPDI_bckgrnd=HPDI(3,:)*1e4; % cases/10,000 people/mnth

% Calculate half-risk distance analytically
if mode_p(2)~=0
    d_half=HalfRskDstnce(p(zthin,:),K0(zthin),typ);
    figure; 
    [mode_d_half,HPDI_d_half]=PlotPstrDistn(d_half,'d_{1/2} (m)',200);
    figure;
    d_half_out=HalfRskDstnce([p(zthin,1:3),zeros(nz,1)],K0(zthin),typ);
    [mode_d_half_out,HPDI_d_half_out]=PlotPstrDistn(d_half_out,'d_{1/2,out} (m)',200);
else
    mode_d_half=NaN;
    HPDI_d_half=NaN(1,2);
    mode_d_half_out=NaN;
    HPDI_d_half_out=NaN(1,2);
end

% Calculate within-HH risk increase (WHHRI)
if mode_p(4)==0
    mode_WHHRI=NaN;
    HPDI_WHHRI=NaN(1,2);
else
    [mode_WHHRI,HPDI_WHHRI]=CalcWthnHHriskIncr(zthin,p,K0);
end

%% INFECTION TIMES
if doPlots
    mode_tE=zeros(nI,1);
    nplot=20;
    %% Cases with both onset and treatment times
    OR1=OR(tI(OR)>maxIP);
    for i=1:nplot
        j=find(I==OR1(i));
%         j=randi(nI,1);
        figure;
        mode_tE(j)=PlotInfctnTimePstrDistn(tEs(j,zthin),tI(OR1(i)),r1,p10,j);
%         mode_tE(j)=PlotInfctnTimePstrDistn(tEs(j,zthin),tI(I(j)),r1,p10,j);
        saveas2(gcf,['E' num2str(j) '_' rslts],savePlots)
        saveas2(gcf,['E' num2str(j) '_' rslts '.eps'],savePlots,'epsc')
    end

    %% Cases without onset or treatment times
    for i=1:min(nplot,nNONR)
        j=find(I==NONR(i));
        figure;
        histogram(tEs(j,zthin),'Normalization','pdf','BinMethod','integers'); hold on
        histogram(tIsNONR(i,zthin),'Normalization','pdf','BinMethod','integers');
        histogram(tRsNONR(i,zthin),'Normalization','pdf','BinMethod','integers')
        set(gca,'FontSize',16);
        xlabel('t (month)','FontSize',16)
        ylabel('Density','FontSize',16)
        h1=legend(['$$E_{' num2str(j) '}$$'],['$$I_{' num2str(j) '}$$'],['$$R_{' num2str(j) '}$$']);
        set(h1,'Interpreter','latex')
        saveas2(gcf,['EIR' num2str(j) '_' rslts],savePlots)
        saveaspdf(gcf,['EIR' num2str(j) '_' rslts])
    end
    
    %% Cases without onset times
    for i=1:nRNO
        j=find(I==RNO(i));
        figure;
        histogram(tEs(j,zthin),'Normalization','pdf','BinMethod','integers'); hold on
        histogram(tIsRNO(i,zthin),'Normalization','pdf','BinMethod','integers');
        set(gca,'FontSize',16);
        xlabel('t (month)','FontSize',16)
        ylabel('Density','FontSize',16)
        h2=legend(['$$E_{' num2str(j) '}$$'],['$$I_{' num2str(j) '}$$']);
        set(h2,'Interpreter','latex')
        saveas2(gcf,['EI' num2str(j) '_' rslts],savePlots)
        saveaspdf(gcf,['EI' num2str(j) '_' rslts])
    end
    %% Cases without treatment times
    mode_tR_ONR=zeros(nONR,1);
    for i=1:nONR
        j=find(I==ONR(i));
        figure; 
        [mode_tE(j),hE]=PlotInfctnTimePstrDistn(tEs(j,zthin),tI(ONR(i)),r1,p10,ONR(i));
        hold on        
        [mode_tR_ONR(i),hR]=PlotRcvryTimePstrDistn(tRsONR(i,zthin),tI(ONR(i)),r0,p0,ONR(i));
        set(gca,'FontSize',16);
        xlabel('t (month)','Interpreter','tex','FontSize',16)
        ylabel('Density','FontSize',16)
        xlim([min(hE.BinEdges) max(hR.BinEdges)])
        h3=legend([hE hR],['$$E_{' num2str(j) '}$$'],['$$R_{' num2str(j) '}$$']);
        set(h3,'Interpreter','latex')
        saveas2(gcf,['ER' num2str(j) '_' rslts],savePlots)
        saveaspdf(gcf,['ER' num2str(j) '_' rslts])
    end
end

%% INCUBATION PERIODS
% Mean IP based on IP distn parameter
mean_IP=mean(r1*(1-p1(zthin))./p1(zthin))+1;
figure;
[mode_IP,HPDI_IP]=PlotPstrDistn(r1*(1-p1(zthin))./p1(zthin)+1,'IP (months)',50);

% Calculate vector of mean IPs for MCMC samples
mean_IPs=mean(IPs,1);
% Plot auto-correlation fn for mean incubation period
figure;
acf(mean_IPs(zthin)',min(200,nz-1));
title('mean IP')

% Plot correlation between mean incubation period and p1
figure; plot(p1(zthin),mean_IPs(zthin),'.')
xlabel('p'); ylabel('mean IP')

%% ASYMPTOMATIC INFECTION AND RECOVERY TIMES
% Plot epi curve with asymptomatic infection
edges=0.5:tmax+0.5;
N1=histcounts(tAs(:,zthin),edges)/nz;
N2=histcounts(tI(I),edges);
N3=histcounts(tP,edges);

t=2002+(0:tmax-1)/12;
clrs=[0.96 0.9 0.8;[254 195 87]/255;[245 150 79]/255;0.8 0.255 0.145;[81 130 187]/255];
figure; plot(t,N1,'Color',clrs(2,:),'LineWidth',1); hold on
plot(t,N2,'Color',clrs(4,:),'LineWidth',1)
plot(t,N3,'Color',clrs(5,:),'LineWidth',1);
set(gca,'FontSize',14)
xlabel('Year'); ylabel('Number')
legend('Asx','VL','PKDL')
saveas(gcf,['EpiCurve_' rslts])
saveas(gcf,['EpiCurve_' rslts '.eps'],'epsc')

% Plot probability of asymptomatic infection before start of study against age
SusAPANIM=setdiff(SusA,IM);
figure; plot(age(SusAPANIM)/12,sum(tAs(SusAPANIM,zthin)==0&tRAs(SusAPANIM,zthin)==0,2)/nz,'.'); hold on
aa=linspace(1,max(age),100)';
plot(aa/12,1-sum(ProbInitStatus(aa,lambda0,p2),2),'LineWidth',2)
set(gca,'FontSize',14)
xlabel('Age (years)'); ylabel('Prob. initially recovered from asymptomatic infection')
legend('imputed prob','model','Location','northwest')
saveas(gcf,['ProbInitRcrvdAsxVsAge_' rslts])
saveas(gcf,['ProbInitRcrvdAsxVsAge_' rslts '.eps'],'epsc')

% Plot posterior probability of being asymptomatically infected during study
figure; set(gcf,'Units','Normalized','OuterPosition',[0 0 1 1]);
bar(sum(tAs(:,zthin)>0 & tAs(:,zthin)<tmax+1,2)/nz); hold on; 
plot(I,ones(nI,1),'r.'); 
ylim([0 1.1])
set(gca,'FontSize',18)
xlabel('Individual'); ylabel('Probability asymptomatically infected during study')
saveas(gcf,['ProbAsxInfctd_' rslts])
saveas(gcf,['ProbAsxInfctd_' rslts '.eps'],'epsc')

% Plots of posterior densities of asymptomatic infection times
for i=[101,7850,20448,22662] % individual 7850 externally immigrated in month 53 into a household with an active KA case
figure; PlotAsxInfctnTimePstrDistn(tAs(i,zthin),probA(i,:),tmax,i)
saveas(gcf,['A' num2str(i) '_' rslts])
saveas(gcf,['A' num2str(i) '_' rslts '.eps'],'epsc')
end

%% PROPOSAL SCALING CONSTANT
figure; plot(0:niters,c,'LineWidth',1);
set(gca,'FontSize',14)
xlabel('Iteration, $$k$$','Interpreter','latex','FontSize',16); ylabel('$$c_k$$','Interpreter','latex','FontSize',16)
saveas(gcf,['SclngFctr_' rslts])
saveas(gcf,['SclngFctr_' rslts '.eps'],'epsc')
