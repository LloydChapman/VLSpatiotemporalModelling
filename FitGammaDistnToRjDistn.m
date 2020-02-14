function [pars,CI]=FitGammaDistnToRjDistn(meanRj)
[pars,CI]=gamfit(meanRj);
figure; histogram(meanRj,'Norm','pdf'); hold on; 
RR=0:0.1:round(max(meanRj)); 
plot(RR,gampdf(RR,pars(1),pars(2)))