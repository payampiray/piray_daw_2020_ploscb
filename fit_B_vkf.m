function [loglik, output] = fit_B_vkf(params,data)
lux      = @(x)(1/(1+exp(-x)));
lambda   = lux(params(1));
v0       = 10*lux(params(2));
omega    = exp(params(3));
params_response = params(4:end);

image    = data.image;
choice   = data.choice;
outcome  = data.outcome;

nt = size(outcome,1);
dv = nan(nt,1);
p  = nan(nt,1);
for i=1:2
    tt = image==i;
    [mu] = vkf_bin(outcome(tt),lambda,v0,omega);
    dv(tt) = mu;
    p(tt) = 1./(1+exp(-mu));
end

Y        = choice==1;
[loglik, beta]   = fit_B_response(p,data.trialvalue,Y,params_response);

if nargout>1
    output = struct('image',image,'outcome',outcome,'lambda',lambda,'v0',v0,'omega',omega,'beta',beta);
end

end
