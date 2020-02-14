function DICminMdl=RunCalcParEstsAndDICdiffs2(rslts,np,burnin,IPD,str,doPlots,savePlots,sprtParas,runName)
% set(0,'DefaultFigureVisible','off')
nMdls=numel(rslts);

%% Calculate parameter estimates and CIs for different models
thin=1; % no thinning currently to match DIC calcs
mode_p=NaN(nMdls,np);
HPDI=NaN(np,2,nMdls);
mode_p1=NaN(nMdls,1);
HPDI1=NaN(nMdls,2);
pcorr=cell(nMdls,1);
pIPtAcorr=cell(nMdls,1);
mode_sptl=NaN(nMdls,1);
HPDI_sptl=NaN(nMdls,2);
mode_bckgrnd=NaN(nMdls,1);
HPDI_bckgrnd=NaN(nMdls,2);
mode_d_half=NaN(nMdls,1);
HPDI_d_half=NaN(nMdls,2);
mode_d_half_out=NaN(nMdls,1);
HPDI_d_half_out=NaN(nMdls,2);
mode_WHHRI=NaN(nMdls,1);
HPDI_WHHRI=NaN(nMdls,2);
mean_IP=NaN(nMdls,1);
mode_IP=NaN(nMdls,1);
HPDI_IP=NaN(nMdls,2);
MdlParEsts=NaN(nMdls,3*(np+1));
for i=1:nMdls
    fprintf([rslts{i} '\n'])
    [mode_p(i,:),HPDI(:,:,i),mode_p1(i),HPDI1(i,:),pcorr{i},pIPtAcorr{i},mode_sptl(i),HPDI_sptl(i,:),mode_bckgrnd(i),HPDI_bckgrnd(i,:),mode_d_half(i),HPDI_d_half(i,:),mode_d_half_out(i),HPDI_d_half_out(i,:),mode_WHHRI(i),HPDI_WHHRI(i,:),mean_IP(i),mode_IP(i),HPDI_IP(i,:)]=ProcessOutput2(rslts{i},burnin,thin,doPlots,savePlots,sprtParas);
    tmp=[];
    for j=1:np
        tmp=[tmp,mode_p(i,j),HPDI(j,:,i)];
    end
    MdlParEsts(i,:)=[tmp,mode_p1(i),HPDI1(i,:)];
%     save(['MdlParEstsFinal' rslts{i}(6:end)])
end

save(['MdlParEstsFinal' IPD runName])

%% Calculate DIC differences from best-fitting model
DICrslts=cellfun(@(x)[str '_' x],rslts,'UniformOutput',false);
DICs=[];
for i=1:numel(DICrslts)
[DICs,DICmin,DICdiffs,RMLs]=CalcDICdiffs(DICrslts{i},DICs);
end
DICminMdl=find(DICs==DICmin);

%% Output parameter estimates and DICs to file
MdlParEstsAndDICs=[MdlParEsts,DICs,DICdiffs];
save(['MdlParEstsAndDICs' IPD sprtParas],'MdlParEstsAndDICs')
% ord=[1:3,5,4,6];
% PrintModesAndCIsToFile(MdlParEstsAndDICs([ord,end-1,end],[1:end-5,end-1]),['MdlParEstsAndDICs' IPD sprtParas '.txt'])

%% Output acceptance rates
acc_rate=NaN(nMdls,11);
for i=1:nMdls
    acc_rate(i,:)=GetAccRates(rslts{i});
end
save(['AccRates' IPD runName],'acc_rate')
