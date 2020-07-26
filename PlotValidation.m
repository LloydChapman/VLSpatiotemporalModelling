function PlotValidation(rslts,ptrue,p1true,runName,Status)

load(rslts)
z=burnin+1:niters;
nz=numel(z);

nbins=50; 
scrnsz=get(0,'ScreenSize');

%% Plot parameter posteriors against actual values used in simulation
figure;
set(gcf, 'Position', [0 0 scrnsz(3) scrnsz(4)/3.8]); %[0 0 1 1/3], 'Units', 'Norm');%
% Start with p to avoid overplotting delta
subplot(1,5,5)
pos = get(gca, 'Position');
pos(1) = 4/5+0.15*1/5;
pos(2) = 0.18;
pos(3) = 0.8*1/5;
pos(4) = 0.8;
set(gca, 'Position', pos)
histogram(p1(z),nbins,'Norm','prob'); hold on
plot(p1true*[1 1],ylim,'r','LineWidth',1)
xlabel('$$p$$','Fontsize',16,'Interpreter','latex')
% ylabel('Density')
for i=1:nu
    j=u(i);
    subplot(1,5,i)
    pos = get(gca, 'Position');
    pos(1) = (i-1)/5+0.15*1/5;
    pos(2) = 0.18;
    pos(3) = 0.8*1/5;
    pos(4) = 0.8;
    set(gca, 'Position', pos)
    histogram(p(z,j),nbins,'Norm','prob'); hold on
    plot(ptrue(j)*[1 1],ylim,'r','LineWidth',1)
    xlabel(['\' pname{j}],'Fontsize',16)
    if i==1
        ylabel('Density')
    end
end
saveas(gcf,['PSTR_DISTNS' runName])
saveas(gcf,['PSTR_DISTNS' runName '.eps'],'epsc')

%% Calculate estimated incidence of asymptomatic infection from MCMC
incA=zeros(tmax+3,nz);
for i=1:nz
    incA(:,z(i))=histcounts(tAs(:,z(i)),'BinM','int'); % counts up to tmax+2
end
median_incA=median(incA(1:end-1,:),2); % exclude counts for tmax+2
Q95=quantile(incA(1:end-1,:),[0.025 0.975],2);

M_incA=zeros(tmax+2,1); % exclude counts for tmax+2 (indexed at tmax+3)
HPDI_incA=zeros(tmax+2,2);
for i=1:tmax+2
    [M_incA(i),HPDI_incA(i,:)]=CalcModeAndHPDI(incA(i,:),50);
end

%% Determine actual incidence of asymptomatic infection in simulation
% tA1=dlmread('tAs_Para1_1.csv');
% tRA1=dlmread('tRAs_Para1_1.csv');
% Status=dlmread('status_mtrx_para1_1.csv');

% Stat0true=dlmread('init_status_para1_1.csv');
Asxtrue=find(simdata.tA<=tmax);
% % actvAtrue=find(Stat0true==2);
% % prevAtrue=find(Stat0true==7 & ~prevK);
% actvAtrue=find(tA1==0 & tRA1>0);
% prevAtrue=find(tA1==0 & tRA1==0);
tAtrue=NaN(n,1);
tAtrue(Asxtrue)=simdata.tA(Asxtrue);
% Remove asx infection times copied to 2nd observations in simulation code
tAtrue(tAtrue<tIM & INTMIG_IN)=NaN;
% Find individuals that remained susceptible till end of study in simulation
Susendtrue=find(all(ismember(Status,[0,1,8]),2) & ~INTMIG_OUT); % only need to check 2nd observations for internal migrators
tAtrue(Susendtrue)=tmax+1;
tRAtrue=NaN(n,1);
tRAtrue(setdiff(Asxtrue,P))=simdata.tR(setdiff(Asxtrue,P));
tRAtrue(PA)=simdata.tDor(PA);
% Remove asx recovery times copied to 2nd observations in simulation code
tRAtrue((tRAtrue<tIM & INTMIG_IN) | (tRAtrue>=tEM & INTMIG_OUT))=NaN;
tRAtrue(Susendtrue)=tmax+1;
% tAtrue(actvAtrue)=0;
% tRAtrue(actvAtrue)=simdata.tR(actvAtrue);
% tAtrue(prevAtrue)=0;
% tRAtrue(prevAtrue)=0;

% Count numbers of new asx infections in each month
incAtrue=histcounts(tAtrue,'BinM','int');
incRAtrue=histcounts(tRAtrue,'BinM','int');

%% Plot inferred numbers of new asymptomatic infections vs actual
tt=0:tmax+1;
tt2=[tt,fliplr(tt)];
Q95R=[Q95(:,1);flipud(Q95(:,2))]';
figure; fill(tt2,Q95R,[0.9 0.9 0.9],'LineStyle','none'); hold on
plot(tt,median_incA,'k-','LineWidth',1)
% HPDR=[HPDI_incA(:,1);HPDI_incA(:,2)]';
% figure; fill(tt2,HPDR,[0.9 0.9 0.9],'LineStyle','none'); hold on
% plot(tt,M_incA,'k-','LineWidth',1)
plot(tt,incAtrue,'k.','MarkerSize',12);
legend('95% CI','median','simulated data')
set(gca,'FontSize',14)
xlabel('Month'); ylabel('No. new asymptomatic infections')
saveas(gcf,['PstrAsxEpiCurveVsSimData' runName])
saveas(gcf,['PstrAsxEpiCurveVsSimData' runName '.eps'],'epsc')
xlim([1 tmax]); ylim([0 max(Q95(2:end-1,2))+1])
% xlim([1 tmax]); ylim([0 max(HPDI_incA(2:end-1,2))+1])
saveas(gcf,['PstrAsxEpiCurveVsSimDataWthnStudy' runName])
saveas(gcf,['PstrAsxEpiCurveVsSimDataWthnStudy' runName '.eps'],'epsc')

%%
% tE1=simdata.tE;
% tI1=simdata.tI;
% tR1=Inf(n,1);
% tR1(~isinf(tI1))=simdata.tR(~isinf(tI1)); % don't think this is sufficient
