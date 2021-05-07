% function [typinfctnsrt,I1srt,infctnsrt,infctrsrt,srcsrt,gen,dist,time,chain]=ReconstructTrnsmssnTree(infctn,infctr,src,typinfctn,I1,dHH,ib,tI,tP,I,P,infctnA)
function [gen,dist,time,chain]=ReconstructTrnsmssnTree(n,nI1,tmax,infctn,infctr,src,I1,dHH,ib,tI,tP,I,P,infctnA)

infctd=(infctn>0 & infctn<tmax+1);
infctn=infctn(infctd);
[~,i]=sort(infctn);

I1=[I1;find(infctnA>0 & infctnA<tmax+1)];
I1srt=I1(i);
infctr=infctr(infctd);
infctrsrt=infctr(i);
% typinfctnsrt=typinfctn(i);

ni=numel(i);
gen1=zeros(ni,1);
dist1=NaN(ni,1);
time1=NaN(ni,1);
chain1=cell(ni,1);
for j=1:ni
    g=1; % should I set g=0 for infections that come from background?
    chain1{j}(1)=I1srt(j);
    I1srttoj=I1srt(1:j);
    infctrj=infctrsrt(j);
    idx=infctrj;
    chain1{j}(2)=idx;
    while ismember(idx,I1srttoj)
        g=g+1;
        idx=infctrsrt(I1srttoj==idx);
        chain1{j}(g+1)=idx;
    end
    gen1(j)=g;
%     if infctrj~=numel(ib)+1 % initial source of infection is not background
%         dist(j)=dHH(ib(infctrj),ib(I1srt(j)));
%         if ismember(idx,I) % index case had presymptomatic infection or VL
%             time(j)=infctn(I1==I1srt(j))-tI(infctrj);
%         elseif ismember(idx,P) % index case had PKDL
%             time(j)=infctn(I1==I1srt(j))-tP(infctrj);
%         else % index case was asymptomatic
%             time(j)=infctn(I1==I1srt(j))-infctnA(infctrj);
%         end
%     end
    if idx~=numel(ib)+1 % initial source of infection is not background
        dist1(j)=dHH(ib(idx),ib(I1srt(j)));
        if ismember(idx,I) % index case had presymptomatic infection or VL
            time1(j)=infctn(I1==I1srt(j))-tI(idx);
        elseif ismember(idx,P) % index case had PKDL
            time1(j)=infctn(I1==I1srt(j))-tP(idx);
        else % index case was asymptomatic
            time1(j)=infctn(I1==I1srt(j))-infctnA(idx);
        end
    end
end
    
[~,j]=sort(i);
gen1=gen1(j);
dist1=dist1(j);
time1=time1(j);
chain1=chain1(j);

gen=NaN(nI1+n,1);
dist=NaN(nI1+n,1);
time=NaN(nI1+n,1);
chain=cell(nI1+n,1);
gen(infctd)=gen1;
dist(infctd)=dist1;
time(infctd)=time1;
chain(infctd)=chain1;