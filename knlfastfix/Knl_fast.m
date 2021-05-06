function [KHH,K0]=Knl_fast(dHH,dHHsqrd,alpha,beta,n,nHH,ib,f,ftransp,I)

KHH=exp(-dHH/alpha); % since we know it's going to be Exp
% KHH(1:nHH+1:end)=0; % N.B. This was incorrect as it set all same household entries to 0
K0=n/sum(sum(f.*ftransp.*KHH)); % proposed fix
% K0=n/sum(sum(f*f'.*KHH)); % 'original' version
% K0=n/sum(K(:)); % calculate normalisation constant
% KHH=K0*KHH(:,ib(I)); % normalise kernel and multiply by transmission rate constant
% KHH=K0*KHH(:,ib); % normalise kernel and multiply by transmission rate constant
KHH=K0*KHH;

%checks on the size of the objects and their classes
%n 25506
%size(KHH) 963x963
%size(K0) 1x1
%size(dHH) 963x1
% size(alpha) 1x1
% size(beta) 1x1
% size(n) 1x1
% size(f) 963x1
% all are doubles below
% class(KHH)
% class(K0)
% class(dHH)
% class(alpha)
% class(beta)
% class(n)
% class(f)
% error("Throw error on first use of Knl_fast()")
% KHH=beta*K0*KHH; % normalise kernel and multiply by transmission rate constant
% g=accumarray([ib(I);nHH],[ones(numel(I),1);0]);
% g=[accumarray(ib(I),1);zeros(nHH-max(ib(I)),1)];
% rate=repelem(KHH,f,g);
% rate=KHH(ib,ib(I)); % expand out to individual-level kernel
