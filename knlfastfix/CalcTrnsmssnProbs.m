function CalcTrnsmssnProbs(rslts,burnin1,thin,runName)

load(rslts)

% Discard burnin
if nargin==1
    z=burnin+1:niters; % if no burnin supplied, use saved default (round(niters/10))
    thin=1;
    runName='';
else
    z=burnin1+1:niters; % if new burnin supplied, use it
end

nbins=50;
zthin=z(1:thin:end);

% Probabilities of getting VL from different sources within the same HH
sptl=p(zthin,1).*K0(zthin)+p(zthin,4);
[p_II,HPDI_II]=TrnsmssnProb(p(zthin,7),sptl,1,nbins);
[p_P1I,HPDI_P1I]=TrnsmssnProb(p(zthin,7),sptl,h1,nbins);
[p_P2I,HPDI_P2I]=TrnsmssnProb(p(zthin,7),sptl,h2,nbins);
[p_P3I,HPDI_P3I]=TrnsmssnProb(p(zthin,7),sptl,h3,nbins);
[p_AI,HPDI_AI]=TrnsmssnProb(p(zthin,7),sptl,h40,nbins);

% Probabilities of getting asymptomatic infection from different sources within the same HH
[p_IA,HPDI_IA]=TrnsmssnProb(1-p(zthin,7),sptl,1,nbins);
[p_P1A,HPDI_P1A]=TrnsmssnProb(1-p(zthin,7),sptl,h1,nbins);
[p_P2A,HPDI_P2A]=TrnsmssnProb(1-p(zthin,7),sptl,h2,nbins);
[p_P3A,HPDI_P3A]=TrnsmssnProb(1-p(zthin,7),sptl,h3,nbins);
[p_AA,HPDI_AA]=TrnsmssnProb(1-p(zthin,7),sptl,h40,nbins);

% Probability of getting VL from background transmission
[p_bckgrndI,HPDI_bckgrndI]=TrnsmssnProb(p(zthin,7),p(zthin,3),1,nbins);

% Probability of getting asymptomatic infection from background transmission
[p_bckgrndA,HPDI_bckgrndA]=TrnsmssnProb(1-p(zthin,7),p(zthin,3),1,nbins);

% Percentage decrease in risk of infection per 100m from a case HH vs 
% directly outside case HH
[RDWD,HPDI_RDWD]=CalcModeAndHPDI(1-exp(-100./p(zthin,2)),nbins);

% Save probabilities
save(['TrnsmssnProbs' runName],'p_II','HPDI_II','p_P1I','HPDI_P1I','p_P2I','HPDI_P2I','p_P3I','HPDI_P3I','p_AI','HPDI_AI','p_IA','HPDI_IA','p_P1A','HPDI_P1A','p_P2A','HPDI_P2A','p_P3A','HPDI_P3A','p_AA','HPDI_AA','p_bckgrndI','HPDI_bckgrndI','p_bckgrndA','HPDI_bckgrndA','RDWD','HPDI_RDWD')


