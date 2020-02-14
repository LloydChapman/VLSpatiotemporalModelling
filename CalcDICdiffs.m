function [DICs,DICmin,DICdiffs,RMLs]=CalcDICdiffs(str,DICs)
load(str)
DICs=[DICs;DIC];
DICmin=min(DICs);
DICdiffs=DICs-DICmin;
RMLs=RltveMdlLklhd(DICs,DICmin);