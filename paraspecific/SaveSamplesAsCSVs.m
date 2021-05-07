function SaveSamplesAsCSVs(rslts,nsmpls,burnin1,runName,varargin)

load(rslts)

if ~exist('z','var')
    z=burnin1+1:niters;
end

if nargin==4
    iters=randperm(numel(z),nsmpls);
else
    iters=varargin{1};
end

ziters=z(iters);
psmpls=p(ziters,:);
p1smpls=p1(ziters);
tEssmpls=tEs(:,ziters);
tAssmpls=tAs(:,ziters);
tRAssmpls=tRAs(:,ziters);
tRsANONRsmpls=tRsANONR(:,ziters);
tRsAONRsmpls=tRsAONR(:,ziters);

dlmwrite(['p' runName '.csv'],psmpls,'precision','%.16f')
dlmwrite(['p1' runName '.csv'],p1smpls,'precision','%.16f')
dlmwrite(['tEs' runName '.csv'],tEssmpls)
dlmwrite(['tAs' runName '.csv'],tAssmpls)
dlmwrite(['tRAs' runName '.csv'],tRAssmpls)
dlmwrite(['tRsANONR' runName '.csv'],tRsANONRsmpls)
dlmwrite(['tRsAONR' runName '.csv'],tRsAONRsmpls)