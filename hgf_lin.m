function [mu2, bad_sigma, mu3, LR] = hgf_lin(y,alpha,nu,kappa,omega)
[N,Q]  = size(y);
y      = [zeros(1,Q); y];
mu3    = nan(N+1,Q);
mu2    = nan(N+1,Q);
sigma2 = nan(N+1,Q);
sigma3 = nan(N+1,Q);

mu2(1,:)    = 0;
sigma2(1,:) = .1;
mu3(1,:)    = 1;
sigma3(1,:) = 1;

LR     = nan(N+1,Q);
v      = nan(N+1,Q);
v(1,:) = exp(omega + kappa*mu3(1,:));

bad_sigma = 0;
for n  = 2:(N+1)       
    expmu3          = exp(omega + kappa*mu3(n-1,:));
    
    delta1          = y(n,:)- mu2(n-1,:);                                  
    LR(n,:)         = (expmu3+sigma2(n-1,:))./(expmu3+sigma2(n-1,:) + alpha);
    mu2(n,:)        = mu2(n-1,:) + LR(n,:).*delta1;                         % Eq 51
    sigma2(n,:)     = LR(n,:)*alpha;                                        % Eq 50
        
    pihat3          = (sigma3(n-1,:) + nu).^-1;                             % Eq 31
    w2              = expmu3 ./ (expmu3 + sigma2(n-1,:) );                  % Eq 32
    r2              = (expmu3 - sigma2(n-1,:))./(expmu3 + sigma2(n-1,:) );  % Eq 33
    delta2          = (sigma2(n,:) +...
        (mu2(n,:)-mu2(n-1,:)).^2)./(sigma2(n-1,:) + expmu3) - 1;            % Eq 34
    
    pi3             = pihat3 + (kappa^2/2).*w2.*(w2+r2.*delta2);            % Eq 29    
    if any(pi3 <= 0)
        bad_sigma = true;
        return;
    end
    sigma3(n,:)     = pi3.^-1;    
    mu3(n,:)        = mu3(n-1,:) + (kappa/2).*sigma3(n,:).*w2.*delta2;      % Eq 30
      
    
    v(n,:)     = exp(omega + kappa*mu3(n,:));
end


vars = [mu2(:) mu3(:) sigma2(:) sigma3(:)];
if any(isnan(vars(:)))
    bad_sigma = true;
end

end
