function RunCalcDIC(rslts,burnin,ss,str)
%% Calculate DICs for models in parallel
nMdls=numel(rslts);
parpool(min(nMdls,12))
parfor i=1:nMdls
    maxNumCompThreads(max(floor(12/nMdls),1))
    CalcDIC(rslts{i},burnin,ss,str);
end