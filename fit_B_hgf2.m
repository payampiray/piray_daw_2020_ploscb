function [loglik] = fit_B_hgf2(params,data)
ux       = @(x)1./(1+exp(-x));
nu       = ux(params(1));
kappa    = ux(params(2));
omega    = -4;
params_response     = (params(3));

image    = data.image;
choice   = data.choice;
outcome  = data.outcome;

nt = size(outcome,1);
p = nan(nt,1);
sig_neg = nan(1,2);
for i=1:2
    tt = image==i;
    [mu1hat,sig_negi] = hgf_bin(outcome(tt), nu, kappa, omega);
    p(tt) = mu1hat;    
    sig_neg(i) = sum(sig_negi);
end

if any(sig_neg)
    loglik = -10^16;
    return;    
end

Y        = choice==1;
loglik   = fit_B_response(p,data.trialvalue,Y,params_response);
end
