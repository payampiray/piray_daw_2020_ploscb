function [loglik] = fit_B_rl(params,data)
ux       = @(x)(1/(1+exp(-x)));
alpha    = ux(params(1));
params_response = (params(2:end));

choice   = data.choice;
outcome  = data.outcome;
image    = data.image;

actions  = 2-choice;
nt = size(outcome,1);
dv = nan(nt,1);
p   = nan(nt,1);
for i=1:2
    tt = image==i;
    [mu] = model_rl(alpha,actions(tt),outcome(tt));
    dv(tt) = 2*mu-1;
    p(tt) = mu;
end

Y        = choice==1;
loglik   = fit_B_response(p,data.trialvalue,Y,params_response);
end

function dv = model_rl(alpha,actions,outcome)
N     = size(outcome,1);
q     = zeros(2,1);
dv    = nan(N,1);
for t = 1:N
    dv(t,:) = q(1) - q(2); % play minus pass
        
    a       = actions(t,:);
    o       = outcome(t,:);
       
    delta   = o - q(a);
    q(a)    = q(a) + alpha.*delta;    
end

end
