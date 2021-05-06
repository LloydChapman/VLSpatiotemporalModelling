function [Rt,Rts]=CalcRt(Rj,RjA,RjI,RjP,onset,tA,onsetI,maxIP,tmax)
Rt=zeros(1,tmax-maxIP);
Rts=zeros(3,tmax-maxIP);
% Rts=zeros(4,tmax);
for j=1:tmax-maxIP
    t=maxIP+j;
    if ismember(t,onset)
        idx=(onset==t);
        Rt(j)=mean(Rj(idx),'omitnan'); % NaNs omitted here to remove NaNs from PKDL cases w/o prior VL (PA), who have NaN in asx part of Rj, as their cntrbtn is counted in PA part
        Rts(1,j)=sum(RjA(tA==t),'omitnan')/sum(idx);
        Rts(2,j)=sum(RjI(onsetI==t),'omitnan')/sum(idx);
        Rts(3,j)=sum(RjP(onsetI==t),'omitnan')/sum(idx); % FOR NOW indexing out of first nIPNIA columns of RjP
%         Rts(4,j)=sum(RjPA(onsetPA==t),'omitnan')/sum(idx);
    end
end