function [loglik, beta] = fit_B_response(p,trialvalue,Y,params)
beta   = exp(params(1));

X      = p.*trialvalue + (1-p)*(-10);
Xb     = bsxfun(@times,X,beta);
f      = 1./(1+exp(-Xb));
p      = f.*Y + (1-f).*(1-Y);
loglik = sum(sum(log(p+eps)));

end