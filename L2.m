function LL2=L2(lambda,pI,tEm)
LL2=sum(log(1-exp(-pI*lambda(tEm)))); % probabilities of being pre-symptomatically infected in given months