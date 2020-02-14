function SaveSamplesAsCSVs(rslts,nsmpls,burnin1,varargin)

load(rslts)

if ~exist('z','var')
    z=burnin1+1:niters;
end

if nargin==3
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

dlmwrite(['p_' rslts(20:end) '.csv'],psmpls,'precision','%.16f')
dlmwrite(['p1_' rslts(20:end) '.csv'],p1smpls,'precision','%.16f')
dlmwrite(['tEs_' rslts(20:end) '.csv'],tEssmpls)
dlmwrite(['tAs_' rslts(20:end) '.csv'],tAssmpls)
dlmwrite(['tRAs_' rslts(20:end) '.csv'],tRAssmpls)
dlmwrite(['tRsANONR_' rslts(20:end) '.csv'],tRsANONRsmpls)
dlmwrite(['tRsAONR_' rslts(20:end) '.csv'],tRsAONRsmpls)