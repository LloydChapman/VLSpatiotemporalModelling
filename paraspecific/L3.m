function LL3=L3(IP,r,p)
LL3=sum(log(nbinpdf(IP-1,r,p))); % probabilities of incubation period durations