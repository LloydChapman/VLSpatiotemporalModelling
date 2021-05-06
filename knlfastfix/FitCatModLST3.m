function [pars,CI]=FitCatModLST3(data,p2)
%FITCATMODLST Fit variable asymptote catalytic model to LST age prevalence data

%% CALCULATE PREVALENCE OF LST POSITIVITY FOR EACH AGE
% Indicator for no current or previous KA
NI=(isnan(data.PREVKA)&data.KA==0);
% NI=(isnan(data.PREVKA)&(data.KA==0|data.FEV_ONS>data.LSTRDT02));
% Make vector of ages at which individuals were tested with LST
a=3:max(data.AGE02(~isnan(data.LST02)&NI));
% Sum number tested at each age
na=zeros(numel(a),1);
for i=1:numel(a)
    na(i)=sum(data.AGE02==a(i)&~isnan(data.LST02)&NI);
end
% Sum number of LST positives at each age
naLpos=zeros(numel(a),1);
for i=1:numel(a)
    naLpos(i)=sum(data.AGE02==a(i)&data.LST02==1&NI);
end

% Remove ages at which no one was tested
iL=(na~=0);
a=a(iL);
naLpos=naLpos(iL);
na=na(iL);
% Calculate prevalence
prevLpos=naLpos./na;

% figure; plot(a,-log(prevLpos),'x');
% figure; bar(a,prevLpos)

%% CATALYTIC MODELS
data1(1,:)=a;
data1(2,:)=a+1;
data1(3,:)=na;
data1(4,:)=naLpos;

% CONST FOI MODEL
% Set initial guess for lambda0
lambda0=0.01;
[pars,CI]=mle(data1(:),'nloglf',@(params,data,cens,freq)nlogLCatMod(params,data,cens,freq,p2),'start',lambda0,'lowerbound',0,'upperbound',Inf);

% Plot prevalence of LST positivity
aa=linspace(1,max(a),1000)';
% p=1-exp(-pars(1)*aa);
prob0=ProbInitStatus(aa,pars,p2);
p=1-sum(prob0,2);
% figure; plot(a,prevLpos,'x',aa,p,'r','LineWidth',1.5);
% legend('data','model');
% xlabel('Age (years)');
% ylabel('Proportion LST+')
% % saveas(gcf,'LSTAgePrev')
% % saveas(gcf,'LSTAgePrev.eps','epsc')

% save('LSTCatModel2Output','pars','CI','a','prevLpos','aa','p')