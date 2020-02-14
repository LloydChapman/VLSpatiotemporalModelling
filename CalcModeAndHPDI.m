function [M,HPDI]=CalcModeAndHPDI(p,nbins,varargin)
if nargin==2
    [probs,edges]=histcounts(p,nbins,'Normalization','prob');
else
    binmthd=varargin{1};
    [probs,edges]=histcounts(p,'Normalization','prob','BinMethod',binmthd);
end
[~,i]=max(probs);
if all(p==p(1))
    M=p(1);
    HPDI=[p(1) p(1)];
else
    M=(edges(i)+edges(i+1))/2;
    sort_probs=sort(probs,'descend');
    HPDIht_idx=find(cumsum(sort_probs)>=0.95,1);
    HPDIht=sort_probs(HPDIht_idx);
    indcs=find(probs>=HPDIht);
    HPDI=[(edges(indcs(1))+edges(indcs(1)+1))/2 (edges(indcs(end))+edges(indcs(end)+1))/2];
end