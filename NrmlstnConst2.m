function K0=NrmlstnConst2(dHH,dHHsqrd,alpha,typ,n,nHH,f)

if strcmp(typ,'Cauchy')
    KHH=1./(1+dHHsqrd/alpha^2);
elseif strcmp(typ,'Exp')
    KHH=exp(-dHH/alpha);
elseif strcmp(typ,'Const')
    KHH=ones(n);
end

KHH(1:nHH+1:end)=0; % set diagonal (same individual) entries to 0
K0=n/sum(sum(f*f'.*KHH)); % calculate normalisation constant