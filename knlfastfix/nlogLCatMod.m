function NLL=nlogLCatMod(lambda0,data,cens,freq,p2)
%NEGLL_VA Negative log-likelihood function for variable asymptote catalytic model

% Reshape data input from vector to matrix
data=[data(1:4:end),data(2:4:end),data(3:4:end),data(4:4:end)]';

n=data(3,:); 
k=data(4,:); 
a=(data(1,:)+data(2,:))/2;
prob0=ProbInitStatus(a',lambda0,p2);
% p=1-exp(-lambda0*a);
p=1-sum(prob0,2)';
LL_a=k.*log(p)+(n-k).*log(1-p);
NLL=-sum(LL_a);