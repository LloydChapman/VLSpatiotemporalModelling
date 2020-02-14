% function [meanBckgrnd,meanASX,meanPRESX,meanKA,meanPKDL]=CalcCntrbtnMgrtnAsx(rslts,nsmpls,burnin1)
function [M,HPDI,iters]=CalcCntrbtnMgrtnAsx(rslts,nsmpls,burnin1)

db='';
rng=[];

load(rslts)
% load('~/Dropbox/Visceral Leishmaniasis/CarynBernData/2010data/data_final')

% Remove first row (initial values) of p if it's not already been removed
% so that same index can be used for p and missing data
if size(p,1)==niters+1
    p=p(2:end,:);
end
% Overwrite asymptomatic infection and recovery times for PKDL cases w/o
% prior VL, if they were not saved, with values from final iteration
if tAs(PA(1))==tmax+2
    tAs(PA,:)=repmat(tA(PA),1,niters);
    tRAs(PA,:)=repmat(tRA(PA),1,niters);
end

% Load data
% load(db)
load('~/Dropbox/Visceral Leishmaniasis/CarynBernData/2010data/data_final2')
% Select data for para
data=data(ismember(data.PARA,para),:);
% Rename longitude and latitude variables
data.Properties.VariableNames{'HHNEWLNG'}='longitude';
data.Properties.VariableNames{'HHNEWLAT'}='latitude';
% Calculate HH distance matrix
dHH=CalcHHDists(data);
dHHsqrd=[];

% Overwrite d0 in old output to make it work with new Knl_fast function
if size(d0,2)>size(d0,1)
    d0=speye(nHH);
end

% Rename hP as hv if hP exists
if exist('hP','var')
    hv=hP;
end

if ~exist('z','var')
    z=burnin1+1:niters;
end

iters=randperm(numel(z),nsmpls);

hP=zeros(nIPNIA,tmax);
for i=1:nIandP
    j=IandP(i);
    hP(I==j,max(tIM(j),tP(j))+1:min(min(tRP(j),tEM(j)),tmax))=hv(j);
end
for i=1:nPI
    j=PI(i);
    hP(nI+i,tP(j)+1:min(min(tRP(j),tEM(j)),tmax))=hv(j);
end
% Need to remove PKDL infectiousness for individual with simultaneous KA
% and PKDL here so we don't double count their contribution, as we treat 
% them as having KA infectiousness until they were treated for KA in MCMC
% code
PpreR=find(tR>tP&~isinf(tR));
hP(ismember(I,PpreR),tP(PpreR)+1:tR(PpreR))=0;
for i=1:nIMP
    j=IMP(i);
    hP(nI+nPI+nA+nIMI+i,tIM(j)+1:min(tRP(j),tmax))=hv(j);
end

% Bckgrnd=zeros(nsmpls,tmax);
% ASX=zeros(nsmpls,tmax);
% PRESX=zeros(nsmpls,tmax);
% KA=zeros(nsmpls,tmax);
% PKDL=zeros(nsmpls,tmax);
nsrcs=5;
FOI=zeros(nsmpls,tmax,nsrcs);

% BckgrndonS=zeros(nsmpls,tmax);
% ASXonS=zeros(nsmpls,tmax);
% PRESXonS=zeros(nsmpls,tmax);
% KAonS=zeros(nsmpls,tmax);
% PKDLonS=zeros(nsmpls,tmax);
FOIonS=zeros(nsmpls,tmax,nsrcs);

nE=sum(sum(tEm));
% Bckgrndm=zeros(nE,nsmpls);
% ASXm=zeros(nE,nsmpls);
% PRESXm=zeros(nE,nsmpls);
% KAm=zeros(nE,nsmpls);
% PKDLm=zeros(nE,nsmpls);
% lambdaEm=zeros(nE,nsmpls);
FOIonE=zeros(nsmpls,nE,nsrcs);

% Bckgrndprop=zeros(nsmpls,1);
% ASXprop=zeros(nsmpls,1);
% PRESXprop=zeros(nsmpls,1);
% KAprop=zeros(nsmpls,1);
% PKDLprop=zeros(nsmpls,1);

% Sus=zeros(nsmpls,tmax);
I2=ismember(I,I1);
% figure(1);
for k=1:nsmpls
    m=z(iters(k));
    
    epsilon=p(m,3);
    
    tE=uint32(tEs(:,m)); % N.B. tEs has dimensions #KA cases x #iterations
    tEm=false(n,tmax);
    tEm((tE(I2)-1)*n+uint32(I1))=1;
    %     for i=1:nI
    %         if tI(I(i))>maxIP
    %             tEm(I(i),tE(i))=1;
    %         end
    %     end
    tA=uint32(tAs(:,m));
    tRA=uint32(tRAs(:,m));
    Asx=uint32(find(tA>=1 & tA<=tmax));
    tAm=false(n,tmax);
    tAm((tA(Asx)-1)*n+Asx)=1;
    prevA=uint32(find(tA==0 & tRA==0));
    actvA=uint32(find(tA==0 & tRA>0));
    IM_INprevAactvA=IM_IN(ismember(IM_OUT,[prevA;actvA]));
    
    if ~inclLST
        S=1-max(preB,preIM)-max(max(max(max(cumsum(tEm,2),cumsum(tAm,2)),cumsum(tPm,2)),cumsum(tDm,2)),cumsum(tEMm,2)); % don't remove LST+ individuals
    else
        S=1-max(preB,preIM)-max(max(max(max(max(cumsum(tEm,2),cumsum(tAm,2)),cumsum(tPm,2)),cumsum(tDm,2)),cumsum(tLm,2)),cumsum(tEMm,2)); % remove LST+ individuals from susceptibles
    end
    S(prevK,:)=0; % remove previous KA cases from susceptibles
    S(prevA,:)=0; % remove previously asymptomatically infected individuals from susceptibles
    S(actvA,:)=0; % remove initially actively asymptomatically infected individuals from susceptibles
    S(tA==0,:)=0; % remove previous and active asymptomatic infections from susceptibles
    S(tI<=maxIP,:)=0; % remove cases with onset before maxIP (also excludes cases with active KA at start of study and a couple of cases with onset in 2002 before immigration)
    S(IpreEXTIM,:)=0; % remove KA cases with onset before or at migration in
    S(EXTIMsoonI,:)=0; % remove KA cases with onset within 6 months of migration in
    S(IpreINTIM,:)=0; % remove KA cases with onset before internal migration in
    S(PpreINTIM,:)=0; % remove PKDL cases with onset before internal migration in 
    S(PpreEXTIM,:)=0; % remove PKDL cases (without prior KA) with onset before external migration in
    S(IM_INprevAactvA,:)=0; % remove susceptibility from 2nd observations of internal migrators asymptomatically infected before the start of the study
    S(IM_IN(ismember(IM_OUT,Asx)),:)=0; % remove susceptibility from 2nd observations of internal migrators asymptomatically infected during 1st observation

    tI(NONR)=tIsNONR(:,m);
    tRorD(NONR)=tRsNONR(:,m);
    tI(RNO)=tIsRNO(:,m);
    tRorD(ONR)=tRsONR(:,m);
    tI(ANONR)=tIsANONR(:,m);
    tRorD(ANONR)=tRsANONR(:,m);
    tRorD(AONR)=tRsAONR(:,m);
    tRL(RLO)=tRLsRLO(:,m);
    tRLR(RLO)=tRLRsRLO(:,m);
    tRL(RLNO)=tRLsRLNO(:,m);
    tRLR(RLNO)=tRLRsRLNO(:,m);
    hE=zeros(nIPNIA,tmax);
    hI=zeros(nIPNIA,tmax);
    for i=1:nI
        j=I(i);
        hE(i,max(max(0,tIM(j)),tE(i))+1:tI(j))=h0;
        hI(i,max(tIM(j),tI(j))+1:min(tEM(j),tRorD(j)))=1;
    end
    for i=1:nA
        hI(nI+nPI+i,1:tRorD(A(i)))=1;
    end
    for i=1:nIMI
        j=IMI(i);
        hI(nI+nPI+nA+i,tIM(j)+1:tRorD(j))=1;
    end
    for i=1:nRL
        j=RL(i);
        hI(IPNIA==j,tRL(j)+1:min(min(tRLR(j),tEM(j)),tmax))=1;
    end
    hA=zeros(n,tmax);
    A1=find(tA>=0 & tA<tmax+1 & tRA>0 & isnan(tP));
    for i=1:numel(A1)
        j=A1(i);
        hA(j,tA(j)+1:min(min(min(tRA(j),tEM(j)),tD(j)),tmax))=p(m,6);
    end
    RAobs1actvA=find(tA==0 & tRA>=tEM & INTMIG_OUT);
    RAobs2actvA=IM_IN(ismember(IM_OUT,RAobs1actvA));    
    RAobs1=find(tA>0 & tA<tEM & tRA>=tEM & INTMIG_OUT);
    RAobs2=IM_IN(ismember(IM_OUT,RAobs1));
    for i=1:numel(RAobs2actvA)
        j=RAobs1actvA(i);
        j1=RAobs2actvA(i);
        hA(j1,rng(j,2)+1:min(min(tRA(j),rng(j1,2)),tmax))=p(m,6);
    end
    for i=1:numel(RAobs2)
        j=RAobs1(i);
        j1=RAobs2(i);
        hA(j1,rng(j,2)+1:min(min(tRA(j),rng(j1,2)),tmax))=p(m,6);
    end
    hA=sparse(hA);
    hPA=zeros(nPA,tmax);
    for i=1:nPA
        j=PA(i);
        hPA(i,max(1,tA(j)+1):tRA(j))=p(m,6);
        hPA(i,tP(j)+1:min(min(tRP(j),tEM(j)),tmax))=hv(j);
    end
    
    KHH=Knl_fast(dHH,dHHsqrd,p(m,2),p(m,1),typ,n,nHH,ib,f,1:n);
    rateHHA=p(m,1)*KHH;
    if ismember(4,u)
        rateHHA=rateHHA+p(m,4)*d0;
    end
    rateHH=rateHHA(:,ib(IPNIA));
    rateHHPA=rateHHA(:,ib(PA));
    lambdaHHA=rateHHA(:,ib)*hA;
    lambdaHHE=rateHH*hE;
    lambdaHHI=rateHH*hI;
    lambdaHHP=rateHH*hP+rateHHPA*hPA;
    lambdaA=lambdaHHA(ib,:);
    lambdaE=lambdaHHE(ib,:);
    lambdaI=lambdaHHI(ib,:);
    lambdaP=lambdaHHP(ib,:);
    lambdaAE=lambdaA(tEm);
    lambdaEE=lambdaE(tEm);
    lambdaIE=lambdaI(tEm);
    lambdaPE=lambdaP(tEm);
% %     lambdaAA=lambdaA(tAm);
% %     lambdaEA=lambdaE(tAm);
% %     lambdaIA=lambdaI(tAm);
% %     lambdaPA=lambdaP(tAm);
%     lambdaEtot=lambdaAE+lambdaIE+lambdaPE+epsilon;
% %     figure(1); %plot(lambdaPE./lambdaE); pause(0.01)
% %     plot(movmean(lambdaPE./lambdaE,10))
    
    % Total FOIs from background, KA, PKDL
%     Bckgrnd(k,:)=n*epsilon; % SHOULD IT BE n OR THE ACTUAL NO. OF PEOPLE (I.E. NOT INCL. UNBORN INDIVIDUALS, DEAD INDIVIDUALS, 2ND OBS FOR INTERNAL MIGRATORS ETC.)
%     ASX(k,:)=sum(lambdaA);
%     PRESX(k,:)=sum(lambdaE);
%     KA(k,:)=sum(lambdaI);
%     PKDL(k,:)=sum(lambdaP);
    FOI(k,:,1)=n*epsilon;
    FOI(k,:,2)=sum(lambdaA);
    FOI(k,:,3)=sum(lambdaE);
    FOI(k,:,4)=sum(lambdaI);    
    FOI(k,:,5)=sum(lambdaP);
    
    % Total FOIs on susceptibles S from background, KA, PKDL
%     BckgrndonS(k,:)=sum(epsilon*S);
%     ASXonS(k,:)=sum(lambdaA.*S);
%     PRESXonS(k,:)=sum(lambdaE.*S);
%     KAonS(k,:)=sum(lambdaI.*S);
%     PKDLonS(k,:)=sum(lambdaP.*S);
    FOIonS(k,:,1)=sum(epsilon*S);
    FOIonS(k,:,2)=sum(lambdaA.*S);
    FOIonS(k,:,3)=sum(lambdaE.*S);
    FOIonS(k,:,4)=sum(lambdaI.*S);
    FOIonS(k,:,5)=sum(lambdaP.*S);
    
    % FOIs on KA cases at their infection times
%     Bckgrndm(:,k)=epsilon;
%     ASXm(:,k)=lambdaAE;
%     PRESXm(:,k)=lambdaEE;
%     KAm(:,k)=lambdaIE;
%     PKDLm(:,k)=lambdaPE;
%     lambdaEm(:,k)=lambdaEtot;
    FOIonE(k,:,1)=epsilon;
    FOIonE(k,:,2)=lambdaAE;
    FOIonE(k,:,3)=lambdaEE;
    FOIonE(k,:,4)=lambdaIE;
    FOIonE(k,:,5)=lambdaPE;
    
%     % Proportion of FOIs on KA cases at infection times from background, KA, PKDL
%     sumlambdaEtot=sum(lambdaEtot);
%     Bckgrndprop(k)=nI*epsilon/sumlambdaEtot;
%     ASXprop(k)=sum(lambdaAE)/sumlambdaEtot;
%     PRESXprop(k)=sum(lambdaEE)/sumlambdaEtot;
%     KAprop(k)=sum(lambdaIE)/sumlambdaEtot;
%     PKDLprop(k)=sum(lambdaPE)/sumlambdaEtot;
    
%     Sus(k,:)=sum(S);
end

%% PLOTS
% Define plot colours
clrs=[0.96 0.9 0.8;[254 195 87]/255;[245 150 79]/255;0.8 0.255 0.145;[81 130 187]/255];
nclrs=size(clrs,1);
% Make time vector (in years)
t=startyr+(0:tmax-1)/12;

% figure; plot(t,mean(Sus))
% xlabel('Time'); ylabel('No. susceptibles')
% saveas(gcf,'NumSscptblesOverTime.png','png')

%% Plot FOI on whole population
nbins=50;
% M_Bckgrnd=zeros(1,tmax);
% M_ASX=zeros(1,tmax);
% M_PRESX=zeros(1,tmax);
% M_KA=zeros(1,tmax);
% M_PKDL=zeros(1,tmax);
% HPDI_Bckgrnd=zeros(tmax,2);
% HPDI_ASX=zeros(tmax,2);
% HPDI_PRESX=zeros(tmax,2);
% HPDI_KA=zeros(tmax,2);
% HPDI_PKDL=zeros(tmax,2);
% for i=1:tmax
%     [M_Bckgrnd(i),HPDI_Bckgrnd(i,:)]=CalcModeAndHPDI(Bckgrnd(:,i),nbins);    
%     [M_ASX(i),HPDI_ASX(i,:)]=CalcModeAndHPDI(ASX(:,i),nbins);
%     [M_PRESX(i),HPDI_PRESX(i,:)]=CalcModeAndHPDI(PRESX(:,i),nbins);
%     [M_KA(i),HPDI_KA(i,:)]=CalcModeAndHPDI(KA(:,i),nbins);
%     [M_PKDL(i),HPDI_PKDL(i,:)]=CalcModeAndHPDI(PKDL(:,i),nbins);
% end
% 
% figure;
% hf=area(t,M);
% % hf=area(t,[mean(Bckgrnd)',mean(ASX)',mean(PRESX)',mean(KA)',mean(PKDL)']);
% for i=1:nclrs
%     hf(i).FaceColor=clrs(i,:);
% end
% xlim([t(1) t(end)])
% set(gca,'FontSize',14)
% xlabel('Time'); ylabel('Total FOI')
% legend('Bckgrnd','Asx','Presx','VL','PKDL','Location','Northwest')
% % saveas(gcf,'AbsltCntrbtnFOIAsxPara1.png','png')
% 
% figure;
% FOItot=ASX+PRESX+KA+PKDL+Bckgrnd;
% hf1=area(t,[mean(Bckgrnd./FOItot)',mean(ASX./FOItot)',mean(PRESX./FOItot)',mean(KA./FOItot)',mean(PKDL./FOItot)']);
% for i=1:nclrs
%     hf1(i).FaceColor=clrs(i,:);
% end
% xlim([t(1) t(end)])
% ylim([0 1])
% set(gca,'FontSize',14)
% xlabel('Time'); ylabel('Relative contribution to FOI')
% legend('Bckgrnd','Asx','Presx','VL','PKDL')
% % saveas(gcf,'RltveCntrbtnFOIAsxPara1.png','png')

PlotCntrbtnModeAndHPDI(FOI,nbins,nsrcs,t,clrs,'Time','FOI (mnth^{-1})',true)
% saveas(gcf,'AbsltCntrbtnFOIAsxPara1HPDI.png','png')
saveas(gcf,'AbsltCntrbtnFOIAsx')
saveas(gcf,'AbsltCntrbtnFOIAsx.png')
PlotCntrbtnModeAndHPDI(bsxfun(@rdivide,FOI,sum(FOI,3)),nbins,nsrcs,t,clrs,'Time','Relative contribution to FOI',true)
% saveas(gcf,'RltveCntrbtnFOIAsxPara1HPDI.png','png')
saveas(gcf,'RltveCntrbtnFOIAsx')
saveas(gcf,'RltveCntrbtnFOIAsx.png')

%% Plot FOI on susceptibles
% figure;
% hf2=area(t,[mean(BckgrndonS)',mean(ASXonS)',mean(PRESXonS)',mean(KAonS)',mean(PKDLonS)']);
% for i=1:nclrs
%     hf2(i).FaceColor=clrs(i,:);
% end
% xlim([t(1) t(end)])
% set(gca,'FontSize',14)
% xlabel('Time'); ylabel('Total FOI')
% legend('Bckgrnd','Asx','Presx','VL','PKDL')
% % saveas(gcf,'AbsltCntrbtnFOIonSAsxPara1.png','png')
%
% figure;
% FOIonS=ASXonS+PRESXonS+KAonS+PKDLonS+BckgrndonS;
% hf3=area(t,[mean(BckgrndonS./FOIonS)',mean(PRESXonS./FOIonS)',mean(ASXonS./FOIonS)',mean(KAonS./FOIonS)',mean(PKDLonS./FOIonS)']);
% for i=1:nclrs
%     hf3(i).FaceColor=clrs(i,:);
% end
% xlim([t(1) t(end)])
% ylim([0 1])
% set(gca,'FontSize',14)
% xlabel('Time'); ylabel('Relative contribution to FOI')
% legend('Bckgrnd','Asx','Presx','VL','PKDL')
% % saveas(gcf,'RltveCntrbtnFOIonSAsxPara1.png','png')

PlotCntrbtnModeAndHPDI(FOIonS,nbins,nsrcs,t,clrs,'Time','FOI x S (mnth^{-1})',true)
% saveas(gcf,'AbsltCntrbtnFOIonSAsxPara1HPDI.png','png')
saveas(gcf,'AbsltCntrbtnFOIonSAsx')
saveas(gcf,'AbsltCntrbtnFOIonSAsx.png')
PlotCntrbtnModeAndHPDI(bsxfun(@rdivide,FOIonS,sum(FOIonS,3)),nbins,nsrcs,t,clrs,'Time','Relative contribution to FOI x S',true)
% saveas(gcf,'RltveCntrbtnFOIonSAsxPara1HPDI.png','png')
saveas(gcf,'RltveCntrbtnFOIonSAsx')
saveas(gcf,'RltveCntrbtnFOIonSAsx.png')

%% Plot FOIs on cases at infection times
% Bckgrndindvdl=mean(Bckgrndm,2);
% ASXindvdl=mean(ASXm,2);
% PRESXindvdl=mean(PRESXm,2);
% KAindvdl=mean(KAm,2);
% PKDLindvdl=mean(PKDLm,2);
% figure; hf4=area([Bckgrndindvdl,PRESXindvdl,ASXindvdl,KAindvdl,PKDLindvdl]);
% for i=1:nclrs
%     hf4(i).FaceColor=clrs(i,:);
% end
% xlim([1 nE])
% set(gca,'FontSize',14)
% xlabel('Case number (by onset)'); ylabel('FOI')
% legend('Bckgrnd','Asx','Presx','VL','PKDL','Location','Northwest')
% % saveas(gcf,'AbsltIndvdlCntrbtnAsxPara1.png','png')
%
% RltveBckgrndindvdl=mean(Bckgrndm./lambdaEm,2);
% RltveASXindvdl=mean(ASXm./lambdaEm,2);
% RltvePRESXindvdl=mean(PRESXm./lambdaEm,2);
% RltveKAindvdl=mean(KAm./lambdaEm,2);
% RltvePKDLindvdl=mean(PKDLm./lambdaEm,2);
% figure; hf5=area([RltveBckgrndindvdl,RltveASXindvdl,RltvePRESXindvdl,RltveKAindvdl,RltvePKDLindvdl]);
% for i=1:nclrs
%     hf5(i).FaceColor=clrs(i,:);
% end
% xlim([1 nE]); ylim([0 1]);
% set(gca,'FontSize',14)
% xlabel('Case number (by onset)'); ylabel('Relative contribution to FOI')
% legend('Bckgrnd','Asx','Presx','VL','PKDL')
% % saveas(gcf,'RltveIndvdlCntrbtnAsxPara1.png','png')

totFOIonE=sum(FOIonE,3);
% figure; area(reshape(mean(FOIonE,1),nE,nsrcs))
% figure; plot(reshape(mean(FOIonE,1),nE,nsrcs))
PlotCntrbtnModeAndHPDI(FOIonE,nbins,nsrcs,1:nE,clrs,'Case number (by onset)','FOI on case (mnth^{-1})',false)
% saveas(gcf,'AbsltIndvdlCntrbtnAsxPara1HPDI.png','png')
saveas(gcf,'AbsltIndvdlCntrbtnAsx')
saveas(gcf,'AbsltIndvdlCntrbtnAsx.png')
PlotCntrbtnModeAndHPDI(bsxfun(@rdivide,FOIonE,totFOIonE),nbins,nsrcs,1:nE,clrs,'Case number (by onset)','Relative contribution to FOI',false)
% saveas(gcf,'RltveIndvdlCntrbtnAsxPara1HPDI.png','png')
saveas(gcf,'RltveIndvdlCntrbtnAsx')
saveas(gcf,'RltveIndvdlCntrbtnAsx.png')

%% Overall contribution to transmission over the course of the epidemic
% figure; hold on
% hf6=histogram(Bckgrndprop);
% hf6.FaceColor=clrs(1,:);
% hf7=histogram(ASXprop);
% hf7.FaceColor=clrs(2,:);
% hf8=histogram(PRESXprop);
% hf8.FaceColor=clrs(3,:);
% hf9=histogram(KAprop);
% hf9.FaceColor=clrs(4,:);
% hf10=histogram(PKDLprop);
% hf10.FaceColor=clrs(5,:);
% legend('Bckgrnd','Asx','Presx','VL','PKDL')
% hold off
% % saveas(gcf,'RltveOvrlCntrbtnAsxPara1.png','png')
% % figure; plot(1:nsmpls,Bckgrnd,1:nsmpls,KA,1:nsmpls,PKDL)
% 
% meanBckgrnd=mean(Bckgrndprop);
% meanPRESX=mean(PRESXprop);
% meanASX=mean(ASXprop);
% meanKA=mean(KAprop);
% meanPKDL=mean(PKDLprop);

figure; hold on
M=zeros(nsrcs,1);
HPDI=zeros(nsrcs,2);
for i=1:nsrcs
    hf6=histogram(sum(FOIonE(:,:,i),2)./sum(totFOIonE,2),nbins,'Normalization','pdf');
    hf6.LineStyle='none';
    hf6.FaceColor=clrs(i,:);
    hf6.FaceAlpha=0.8;
    [M(i),HPDI(i,:)]=CalcModeAndHPDI(sum(FOIonE(:,:,i),2)./sum(totFOIonE,2),nbins);
end
set(gca,'FontSize',14)
xlabel('Relative overall contribution')
ylabel('Density')
legend('Bckgrnd','Asx','Presx','VL','PKDL')
hold off
saveas(gcf,'RltveOvrlCntrbtnAsx')
saveas(gcf,'RltveOvrlCntrbtnAsx.png')

%% Save mode and HPDI for relative overall contribution
save('ModeAndHPDIRltveOvrlCntrbtnAsx','M','HPDI','iters')
