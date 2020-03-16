function [loglik, output] = fit_A_vkf(params,data)

lux      = @(x)(1/(1+exp(-x)));
lambda   = lux(params(1));
v0       = 10*lux(params(2));
omega    = exp(params(3));
params_response = params(4:7);

choice   = data.choice;
outcome  = data.outcome;

outcome(:,1:2) = 2*outcome(:,1:2) - 1;
outcome(:,3:4) = 2*outcome(:,3:4) + 1;
outcome(choice==2) = -outcome(choice==2);

outcome(outcome==-1) = 0;

dv  = vkf_bin(outcome,lambda,v0,omega);


Y        = choice==1;
[loglik, beta]   = fit_A_response(dv,Y,params_response);

if nargout>1
    output = struct('outcome',outcome,'lambda',lambda,'v0',v0,'omega',omega,'beta',beta);
end

end

