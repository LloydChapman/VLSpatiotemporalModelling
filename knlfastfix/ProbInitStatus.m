% function probs=ProbInitStatus(a,lambda,pI,p2)
function probs=ProbInitStatus(a,lambda,p2)
% probs=[exp(-lambda*a),(1-exp(-(1-pI)*lambda))*((1-p2).^a-exp(-lambda*a))/(1-p2-exp(-lambda)),(1-exp(-(1-pI)*lambda))*((1-exp(-lambda*a))/(1-exp(-lambda))-((1-p2).^a-exp(-lambda*a))/(1-p2-exp(-lambda)))];
probs=[exp(-lambda.*a),(1-exp(-lambda)).*((1-p2).^a-exp(-lambda.*a))./(1-p2-exp(-lambda))];
