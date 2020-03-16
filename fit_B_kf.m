function [loglik] = fit_B_kf(params,data)
sigma    = exp(params(1));
omega    = exp(params(2));
params_response = (params(3:end));

image    = data.image;
choice   = data.choice;
outcome  = data.outcome;

nt = size(outcome,1);
dv = nan(nt,1);
p = nan(nt,1);
for i=1:2
    tt = image==i;
    [mu] = kf(outcome(tt),sigma,omega);
    dv(tt) = 2*mu-1;
    p(tt) = mu;
end

Y        = choice==1;
loglik   = fit_B_response(p,data.trialvalue,Y,params_response);
end

function dv = kf(y,sigma,omega)

[nt,nq] = size(y);

m       = .5*ones(1,nq);
w       = sigma*ones(1,nq);

dv      = nan(nt,nq);

for t  = 1:nt      
    dv(t,:)     = (m(1,:));
        
    k           = (w+sigma)./(w+sigma + omega);
    delta       = y(t,:) - m;
    m           = m + k.*delta;
    w           = (1-k).*(w+sigma);
end

end
