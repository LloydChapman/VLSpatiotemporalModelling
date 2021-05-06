function [dists,onsetinfctr,rcvryinfctr,times,SI]=GetDistsAndTimes(infctn,infctr,src,dHH,ib,I1,tA,tI,tP,tRorD,tRP)
% Find distances between most likely infectors and KA cases
nonbckgrnd=(infctr~=numel(ib)+1);
nI1=numel(I1);
% if nI1==numel(nonbckgrnd)
    idx=I1(nonbckgrnd);
% else
%     idx=intersect(I1,find(nonbckgrnd));
% end
dists=NaN(nI1,1);
dists(nonbckgrnd)=dHH(sub2ind(size(dHH),ib(idx),ib(infctr(nonbckgrnd))));

% Find onsets of most likely infectors
onsetinfctr=NaN(nI1,1);
onsetinfctr(nonbckgrnd&src==2)=tA(infctr(nonbckgrnd&src==2));
onsetinfctr(nonbckgrnd&src==3)=tI(infctr(nonbckgrnd&src==3)); % N.B. Will give -ve time difference between onset and secondary infection for cases infected by pre-symptomatics - think abut whether this makes sense
onsetinfctr(nonbckgrnd&src==4)=tI(infctr(nonbckgrnd&src==4));
onsetinfctr(nonbckgrnd&src==6)=tP(infctr(nonbckgrnd&src==6));
% Find treatment times of most likely infectors
rcvryinfctr=NaN(numel(I1),1);
rcvryinfctr(nonbckgrnd&src==3)=tRorD(infctr(nonbckgrnd&src==3));
rcvryinfctr(nonbckgrnd&src==4)=tRorD(infctr(nonbckgrnd&src==4));
rcvryinfctr(nonbckgrnd&src==6)=tRP(infctr(nonbckgrnd&src==6));
% Find times between onsets of most likely infectors and infections of
% KA cases
times=NaN(nI1,1);
times(nonbckgrnd)=infctn(nonbckgrnd)-onsetinfctr(nonbckgrnd);
% Find serial intervals
SI=NaN(nI1,1);
SI(nonbckgrnd)=min(double(tA(idx)),tI(idx))-onsetinfctr(nonbckgrnd);
%     % Find generation times - not so easy as some infectors won't have infection times
%     gt(nonbckgrnd)=infctn(nonbckgrnd)-tE(infctr(nonbckgrnd));
