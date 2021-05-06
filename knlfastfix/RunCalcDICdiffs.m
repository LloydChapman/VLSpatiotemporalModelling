function [DICs,DICmin,DICdiffs,RMLs,DICminMdl]=RunCalcDICdiffs(rslts,str)
DICrslts=cellfun(@(x)[str '_' x],rslts,'UniformOutput',false);
DICs=[];
for i=1:numel(DICrslts)
[DICs,DICmin,DICdiffs,RMLs]=CalcDICdiffs(DICrslts{i},DICs);
end
DICminMdl=find(DICs==DICmin);