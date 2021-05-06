function LL4=L4(lambda,pI,tAm)
LL4=sum(log(1-exp(-(1-pI)*lambda(tAm)))); % probabilities of being asymptomatically infected in given months