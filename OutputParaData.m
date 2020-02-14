clear
set(0,'defaultAxesFontSize',14)
set(0,'defaultLineLineWidth',1)

%% Load data
load('data_final2.mat')

%% Set parameters for calculating incidence
% Set start year and start month
startyr=2002;
startmo=1;
 
% Set end year and month
endyr=2010;
endmo=12;

% Find unique paras in data
para=unique(data.PARA);
npara=numel(para);

% Make STATA month origin
origin=stata_month(startyr,startmo)-1;

% Make indicator variables for 1st and 2nd observations for individuals who internally migrated 
INTMIG_OUT=(data.INTMIG_OUT==1);
INTMIG_IN=(data.INTMIG_IN==1);

% Create logical index vector for duplicate observations of KA for KA cases 
% who migrated internally
KothrObs=((data.KA>=data.MIG_OUT&INTMIG_OUT)|(data.KA<data.MIG_IN&INTMIG_IN));

% Create logical index vector for duplicate observations of PKDL for PKDL 
% cases who migrated internally
PothrObs=((data.PKDL>=data.MIG_OUT&INTMIG_OUT)|(data.PKDL<data.MIG_IN&INTMIG_IN));

% Create indicator variable for KA onset between Jan 2002 and Dec 2010
K02_10=(data.KA_1210&data.KAYR>=startyr&~KothrObs);

% Create indicator variable for PKDL onset between Jan 2002 and Dec 2010
P02_10=(data.PKDL_1210&~PothrObs);

% End time
tmax=stata_month(endyr,endmo)-origin;

% Time vectors
t=1:tmax;
mos=startyr+(t-1)/12;
yrs=startyr:endyr;
nyrs=numel(yrs);

%% Total KA and PKDL cases in each para
tbl=tabulate(data.PARA_ID(data.INTMIG_IN~=1));
paradata=cell2table(tbl(:,1:2),'VariableNames',{'PARA_ID','n'});

% % Count KA cases in each para
% tbl1=tabulate(data.PARA_ID(K02_10));
% paradata.VLcases=cell2mat(tbl1(:,2));
% 
% % Count PKDL cases in each para
% tbl2=tabulate(data.PARA_ID(P02_10));
% paradata=outerjoin(paradata,cell2table(tbl2(:,1:2),'VariableNames',{'PARA_ID','PKDLcases'}),'MergeKeys',true);
% paradata.PKDLcases(isnan(paradata.PKDLcases))=0;

%% Calculate KA and PKDL incidence
% Monthly para-level populations, KA cases and PKDL cases
pop=zeros(npara,tmax);
l=cell(npara,1);
cases=zeros(npara,tmax);
casesP=zeros(npara,tmax);
for i=1:npara
    for j=1:tmax
        mo=origin+j;
        pop(i,j)=sum(data.PARA==para(i) & data.ENTRY<=mo & (data.EXIT>mo | data.EXIT==origin+tmax)); %&(~data.KA_1210|data.KA>=mo)&(~data.PKDL_1210|data.PKDL>=mo)
        cases(i,j)=sum(data.PARA==para(i) & data.KA==mo & K02_10);
        casesP(i,j)=sum(data.PARA==para(i) & data.PKDL==mo & P02_10);
    end
    l{i}=num2str(para(i));
end
figure; plot(mos,pop)
xlabel('Time'); ylabel('Population')
hl=legend(l); title(hl,'Para')
ovrl_pop=sum(pop,1);
figure; plot(mos,ovrl_pop)
xlabel('Year'); ylabel('Total population')

% Monthly para-level KA incidence
inc=1e4*12*cases./pop;
figure; plot(mos,inc)
hl=legend(l); title(hl,'Para')
xlabel('Time'); ylabel('VL incidence (cases/10,000/yr)')
saveas(gcf,'ParaMonthlyKAinc')
saveas(gcf,'ParaMonthlyKAinc.eps','epsc')

% Monthly para-level PKDL incidence
incP=1e4*12*casesP./pop;
figure; plot(mos,incP)
hl=legend(l); title(hl,'Para')
xlabel('Time'); ylabel('PKDL incidence (cases/10,000/yr)')
saveas(gcf,'ParaMonthlyPKDLinc')
saveas(gcf,'ParaMonthlyPKDLinc.eps','epsc')

% Monthly overall KA and PKDL incidence
ovrl_cases=sum(cases,1);
ovrl_casesP=sum(casesP,1);
ovrl_inc=1e4*12*ovrl_cases./ovrl_pop;
ovrl_incP=1e4*12*ovrl_casesP./ovrl_pop;
figure; h=plot(mos,ovrl_inc,mos,ovrl_incP);
h(1).Color=[0.8 0.255 0.145];
h(2).Color=[81 130 187]/255;
set(gca,'FontSize',24)
xlabel('Time'); ylabel('Incidence (cases/10,000/yr)'); legend('VL','PKDL')
saveas(gcf,'MonthlyKAandPKDLinc')
saveas(gcf,'MonthlyKAandPKDLinc.eps','epsc')

% Annual para-level populations, KA cases and PKDL cases
% pop_yr=zeros(npara,nyrs);
cases_yr=zeros(npara,nyrs);
casesP_yr=zeros(npara,nyrs);
for i=1:npara
    for j=1:nyrs
        yr=startyr+j-1;
%         pop_yr(i,j)=sum(data.PARA==para(i) & data.ENTRY_Y<=yr & (data.EXIT_Y>yr | data.EXIT_Y==endyr));
        cases_yr(i,j)=sum(data.PARA==para(i) & data.KAYR==yr & K02_10);
        casesP_yr(i,j)=sum(data.PARA==para(i) & data.PKDLYR==yr & P02_10);
    end
end
% Calculate mean annual populations from pop to account for births, deaths
% and migration
pop_yr=squeeze(mean(reshape(pop,npara,12,nyrs),2));
figure; plot(yrs,pop_yr)
xlabel('Year'); ylabel('Population')
hl=legend(l); title(hl,'Para')
tot_pop_yr=sum(pop_yr,1);
figure; plot(yrs,tot_pop_yr)
xlabel('Year'); ylabel('Total population')

% Annual para-level KA incidence
inc_yr=1e4*cases_yr./pop_yr;
figure; set(axes,'LineStyleOrder',{'-','--','-.',':'}')
hold all
plot(yrs,inc_yr,'LineWidth',1)
hl=legend(l,'FontSize',9); title(hl,'Para')
xlabel('Year'); ylabel('VL incidence (cases/10,000/yr)')
saveas(gcf,'ParaAnnlKAinc')
saveas(gcf,'ParaAnnlKAinc.eps','epsc')

% Annual para-level PKDL incidence
incP_yr=1e4*casesP_yr./pop_yr;
figure; set(axes,'LineStyleOrder',{'-','--','-.',':'}')
hold all
plot(yrs,incP_yr)
hl=legend(l,'FontSize',9,'Location','northwest'); title(hl,'Para')
xlabel('Year'); ylabel('PKDL incidence (cases/10,000/yr)')
saveas(gcf,'ParaAnnlPKDLinc')
saveas(gcf,'ParaAnnlPKDLinc.eps','epsc')

% Annual para-level KA incidence distribution
figure; histogram(inc_yr,20)
xlim([0 Inf])
set(gca,'FontSize',24)
xlabel('VL incidence (cases/10,000/yr)'); ylabel('Number of paras')
saveas(gcf,'ParaKAincDistn')
saveas(gcf,'ParaKAincDistn.eps','epsc')
figure; histogram(incP_yr,20)
xlim([0 Inf])
set(gca,'FontSize',24)
xlabel('PKDL incidence (cases/10,000/yr)'); ylabel('Number of paras')
saveas(gcf,'ParaPKDLincDistn')
saveas(gcf,'ParaPKDLincDistn.eps','epsc')

% Calculate para-level average KA and PKDL incidences
mean_pop=mean(pop,2);
paradata.MeanPop=mean_pop;
tot_cases_yr=sum(cases_yr,2);
tot_casesP_yr=sum(casesP_yr,2);
paradata.VLcases=tot_cases_yr;
paradata.PKDLcases=tot_casesP_yr;
paradata.MeanVLinc=1e4*tot_cases_yr./(mean_pop*nyrs);
paradata.MeanPKDLinc=1e4*tot_casesP_yr./(mean_pop*nyrs);
paradata.MinAnnlVLinc=min(inc_yr,[],2);
paradata.MaxAnnlVLinc=max(inc_yr,[],2);
paradata.MinAnnlPKDLinc=min(incP_yr,[],2);
paradata.MaxAnnlPKDLinc=max(incP_yr,[],2);

%% Determine which cluster paras are in
data.CLUSTER=cell(size(data,1),1);
data.CLUSTER(data.HHNEWLNG<90.275)={'NW'};
data.CLUSTER(data.HHNEWLNG>90.275)={'SE'};
clstr=unique(data(:,{'PARA_ID','CLUSTER'}),'rows');
paradata.CLUSTER=clstr{:,2};
paradata.VILLAGE=cellfun(@(x)x(1:2),paradata.PARA_ID,'UniformOutput',false);
paradata.PARA_NAME=cellfun(@(x)x(3:4),paradata.PARA_ID,'UniformOutput',false);
paradata.PARA=unique(data.PARA);

%% Calculate cumulative incidence
paradata.CumVLinc=paradata.VLcases./paradata.n;
paradata.CumPKDLinc=paradata.PKDLcases./paradata.n;

%% Overall values
n=sum(paradata.n);
tot_ovrl_cases=sum(paradata.VLcases);
tot_ovrl_casesP=sum(paradata.PKDLcases);
mean_ovrl_pop=mean(ovrl_pop);
mean_ovrl_inc=1e4*tot_ovrl_cases/(mean_ovrl_pop*nyrs);
mean_ovrl_incP=1e4*tot_ovrl_casesP/(mean_ovrl_pop*nyrs);
ovrl_pop_yr=mean(reshape(ovrl_pop,12,nyrs),1);
ovrl_cases_yr=sum(cases_yr,1);
ovrl_casesP_yr=sum(casesP_yr,1);
ovrl_inc_yr=1e4*ovrl_cases_yr./ovrl_pop_yr;
ovrl_incP_yr=1e4*ovrl_casesP_yr./ovrl_pop_yr;
paradata(npara+1,{'PARA','PARA_ID','n','MeanPop','VLcases','PKDLcases','MeanVLinc','MeanPKDLinc','MinAnnlVLinc','MaxAnnlVLinc','MinAnnlPKDLinc','MaxAnnlPKDLinc','CLUSTER','VILLAGE','PARA_NAME','CumVLinc','CumPKDLinc'})={NaN,'Total',n,mean_ovrl_pop,tot_ovrl_cases,tot_ovrl_casesP,mean_ovrl_inc,mean_ovrl_incP,min(ovrl_inc_yr),max(ovrl_inc_yr),min(ovrl_incP_yr),max(ovrl_incP_yr),'','','',tot_ovrl_cases/n,tot_ovrl_casesP/n};

%% Write output
paradata=paradata(:,{'PARA','PARA_ID','CLUSTER','VILLAGE','PARA_NAME','n','MeanPop','VLcases','PKDLcases','CumVLinc','CumPKDLinc','MeanVLinc','MeanPKDLinc','MinAnnlVLinc','MaxAnnlVLinc','MinAnnlPKDLinc','MaxAnnlPKDLinc'});
writetable(paradata,'ParaData2.csv');

%% Plot map of paras
npara=max(data.PARA);
figure; hold on
l=cell(npara,1);
cmap=jet(npara);
for i=1:npara
    plot(data.HHNEWLNG(data.PARA==para(i)),data.HHNEWLAT(data.PARA==para(i)),'.','Color',cmap(i,:));
    l{i}=num2str(para(i));
end
legend(l)
xlabel('lon'); ylabel('lat')
hold off
saveas(gcf,'ParaMap')
saveas(gcf,'ParaMap.png')

%% Plot map of villages
data.VILLAGE_ID=cellfun(@(x)x(1:2),data.PARA_ID,'UniformOutput',false);
[~,~,data.VILLAGE]=unique(data.VILLAGE_ID);
nvillage=numel(unique(data.VILLAGE));
figure; hold on
l=cell(nvillage,1);
cmap=jet(nvillage);
for i=1:nvillage
    plot(data.HHNEWLNG(data.VILLAGE==i),data.HHNEWLAT(data.VILLAGE==i),'.','Color',cmap(i,:));
    l{i}=num2str(i);
end
legend(l)
xlabel('lon'); ylabel('lat')
hold off
saveas(gcf,'VillageMap')
saveas(gcf,'VillageMap.png')