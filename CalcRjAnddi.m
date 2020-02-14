function [Rj,di,ti,pij_ti]=CalcRjAnddi(lambdaij,lambdai,d,tij)
% pij=bsxfun(@rdivide,lambdaij,lambdai);
% pij=1-exp(-lambdaij);
% Probability that i was infected by j conditional on i being infected at
% time ti
% pij_ti=bsxfun(@rdivide,pij,pi);
pij_ti=bsxfun(@rdivide,lambdaij,lambdai);
Rj=sum(pij_ti,1);
% pij_ti1=bsxfun(@rdivide,pij,sum(pij,2,'omitnan'));
pij_ti1=bsxfun(@rdivide,lambdaij,sum(lambdaij,2,'omitnan'));
di=sum(d.*pij_ti1,2,'omitnan');
allzero=all(pij_ti1==0,2);
di(allzero)=NaN;
ti=sum(tij.*pij_ti1,2,'omitnan');
ti(allzero)=NaN;