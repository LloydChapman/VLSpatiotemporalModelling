function [KHH,K0]=Knl_fast(dHH,dHHsqrd,alpha,beta,typ,n,nHH,ib,f,I)

if strcmp(typ,'Cauchy')
    KHH=1./(1+dHHsqrd/alpha^2);
elseif strcmp(typ,'Exp')
    KHH=exp(-dHH/alpha);
elseif strcmp(typ,'Const')
    KHH=ones(nHH);
end

% KHH(1:nHH+1:end)=0; % N.B. This was incorrect as it set all same household entries to 0
K0=n/sum(sum(f*f'.*KHH)); % calculate normalisation constant
% K0=n/sum(K(:)); % calculate normalisation constant
% KHH=K0*KHH(:,ib(I)); % normalise kernel and multiply by transmission rate constant
% KHH=K0*KHH(:,ib); % normalise kernel and multiply by transmission rate constant
KHH=K0*KHH;
% KHH=beta*K0*KHH; % normalise kernel and multiply by transmission rate constant
% g=accumarray([ib(I);nHH],[ones(numel(I),1);0]);
% g=[accumarray(ib(I),1);zeros(nHH-max(ib(I)),1)];
% rate=repelem(KHH,f,g);
% rate=KHH(ib,ib(I)); % expand out to individual-level kernel

