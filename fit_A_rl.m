function [loglik] = fit_A_rl(params,data)

ux       = @(x)(1/(1+exp(-x)));
alpha    = ux(params(1));
params_response = params(2:5);

choice   = data.choice;
outcome  = data.outcome;

[dv]     = model_rl(alpha,choice,outcome);

Y        = choice==1;
[loglik] = fit_A_response(dv,Y,params_response);
end

function dv = model_rl(alpha,actions,outcome)
Q      = size(outcome,2);
N      = size(outcome,1);

q       = zeros(2,Q);
dv      = nan(N,Q);

for t = 1:N
    dv(t,:) = q(1,:) - q(2,:);
    
    a      = actions(t,:);    
    o      = outcome(t,:);
    idx    = sub2ind([2 Q],a,(1:Q));
       
    delta  = o - q(idx);
    q(idx) = q(idx) + alpha.*delta;    
end

end
