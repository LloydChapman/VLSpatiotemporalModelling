function RunCalcDIC(rslts,burnin,ss,str)
%% Calculate DICs for models in parallel
nMdls=numel(rslts);
numcores=feature('numcores');
parpool(min(nMdls,numcores))
parfor i=1:nMdls
    maxNumCompThreads(max(floor(numcores/nMdls),1))
    CalcDIC(rslts{i},burnin,ss,str);
end