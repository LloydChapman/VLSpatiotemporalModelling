function DIC=dic(D,MD)
%DIC Deviance information criterion
% DIC=2*mean(D(p))-D(phat), where D(p)=-2*log(L(p)) is the deviance for
% the parameter vector p and likelihood L(p) and phat is a point estimate
% of p, e.g. the mode or mean, from its posterior distribution
DIC=2*MD-D; % MD is the mean deviance from the MCMC and D is the deviance evaluated at the posterior modes