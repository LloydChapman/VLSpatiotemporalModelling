function h=PlotDevDistn(str,burnin1)
% Load MCMC output
load(str,'LL','niters')

% Redefine burnin
z=burnin1+1:niters;

% Plot deviance distribution
h=histogram(-2*LL(z),50,'Normalization','pdf');
xlabel('Posterior deviance'); ylabel('Density')
