function [loglik] = fit_A_kf(params,data)
sigma    = exp(params(1));
omega    = exp(params(2));
params_response = params(3:6);

choice   = data.choice;
outcome  = data.outcome;

outcome(:,1:2) = 2*outcome(:,1:2) - 1;
outcome(:,3:4) = 2*outcome(:,3:4) + 1;
outcome(choice==2) = -outcome(choice==2);

outcome(outcome==-1) = 0;

[dv]  = kf(outcome,sigma,omega);
dv    = 2*dv-1;

Y        = choice==1;
loglik   = fit_A_response(dv,Y,params_response);
end
% 
function dv = kf(y,sigma,omega)

[nt,nq] = size(y);

m       = zeros(1,nq);
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
