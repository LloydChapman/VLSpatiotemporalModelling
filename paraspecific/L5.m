% function LL5=L5(S,age,Sus,actvA,prevA,Pres0,lambda0,pI,p2)
% ageS=age(intersect(union(Sus,find(S(:,1))),Pres0));
function LL5=L5(age,S0,actvA,prevA,lambda0,pI,p2)
ageS=age(S0);
ageS=ageS(ageS~=0); % N.B. need to think about what to do with 83 individuals without DOB; just exclude them from LL5 for now
ageA=age(actvA);
ageA=ageA(ageA~=0); % N.B. need to think about what to do with 83 individuals without DOB; just exclude them from LL6 for now
ageRA=age(prevA);
ageRA=ageRA(ageRA~=0); % N.B. need to think about what to do with 83 individuals without DOB; just exclude them from LL6 for now
% probS=ProbInitStatus(ageS,lambda0,pI,p2);
% probA=ProbInitStatus(ageA,lambda0,pI,p2);
% probRA=ProbInitStatus(ageRA,lambda0,pI,p2);
% % LL5=sum(log(probS(:,1)))+sum(log(probA(:,2)))+sum(log(probRA(:,3))); % probabilities of being initially susceptible, initially asymptomatically infeceted, and previously asymptomatically infected 
% if ~isempty(probA)
%     LL5=sum(log(probS(:,1)))+sum(log(probA(:,2)))+sum(log(1-sum(probRA,2))); % probabilities of being initially susceptible, initially asymptomatically infected, and previously asymptomatically infected 
% else
%     LL5=sum(log(probS(:,1)))+sum(log(1-sum(probRA,2))); % probabilities of being initially susceptible and previously asymptomatically infected 
% end
LL5=-sum(lambda0*ageS)+sum(log((1-exp(-lambda0))*((1-p2).^ageA-exp(-lambda0*ageA))/(1-p2-exp(-lambda0))))+sum(log(1-exp(-lambda0*ageRA)-(1-exp(-lambda0))*((1-p2).^ageRA-exp(-lambda0*ageRA))/(1-p2-exp(-lambda0)))); % probabilities of being initially susceptible, initially asymptomatically infected, and previously asymptomatically infected 