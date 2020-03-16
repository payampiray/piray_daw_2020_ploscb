function [loglik] = fit_A_hgf(params,data)
ux       = @(x)1./(1+exp(-x));
nu       = ux(params(1));
kappa    = ux(params(2));
omega    = params(3) -1;
params_response = (params(4:7));

choice   = data.choice;
outcome  = data.outcome;

outcome(:,1:2) = 2*outcome(:,1:2) - 1;
outcome(:,3:4) = 2*outcome(:,3:4) + 1;
outcome(choice==2) = -outcome(choice==2);

outcome(outcome==-1) = 0;

[mu1hat, sig_neg] = hgf_bin(outcome, nu, kappa, omega);
dv = 2*mu1hat - 1;
if sig_neg    
    loglik = -10^16;
    return;
    dv = zeros(size(dv));
    params_response = 0*params_response;
end

Y        = choice==1;
[loglik]   = fit_A_response(dv,Y,params_response);
end
