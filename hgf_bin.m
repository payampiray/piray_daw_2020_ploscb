function [m, bad_traj, mu2, mu3, sigma2] = hgf_bin(y,nu,kappa,omega)

[N,Q]  = size(y);
y      = [zeros(1,Q); y];
mu3    = nan(N+1,Q);
mu2    = nan(N+1,Q);
sigma2 = nan(N+1,Q);
sigma3 = nan(N+1,Q);
mu1hat = nan(N+1,Q);

mu2(1,:)    = 0;
sigma2(1,:) = 0.1;
mu3(1,:)    = 1;
sigma3(1,:) = 1;

bad_traj = false;

ux = @(x)(1./(1+exp(-x)));
for n  = 2:(N+1)       
    expmu3          = exp(omega + kappa*mu3(n-1,:));
    
    mu1hat(n,:)     = ux(mu2(n-1,:));                                       % Eq 24
    delta1          = y(n,:)- mu1hat(n,:);                                  % Eq 25 
    sigma1hat       = mu1hat(n,:).*(1-mu1hat(n,:));                         % Eq 26
    sigma2hat       = sigma2(n-1,:)+ expmu3;                                % Eq 27
    pi2             = sigma2hat.^-1 + sigma1hat;                            % Eq 28
    
    LR              = pi2.^-1;
    mu2(n,:)        = mu2(n-1,:) + LR.*delta1;                              % Eq 23
    sigma2(n,:)     = LR;                                                   % Eq 22
        
    pihat3          = (sigma3(n-1,:) + nu).^-1;                             % Eq 31
    w2              = expmu3 ./ (expmu3 + sigma2(n-1,:) );                  % Eq 32
    r2              = (expmu3 - sigma2(n-1,:))./(expmu3 + sigma2(n-1,:) );  % Eq 33
    delta2          = (sigma2(n,:) +...
        (mu2(n,:)-mu2(n-1,:)).^2)./(sigma2(n-1,:) + expmu3) - 1;            % Eq 34
    
    pi3             = pihat3 + (kappa^2/2).*w2.*(w2+r2.*delta2);            % Eq 29
    if any(pi3 <= 0)
        bad_traj = true;
    end
    sigma3(n,:)     = pi3.^-1;        
    mu3(n,:)        = mu3(n-1,:) + (kappa/2).*sigma3(n,:).*w2.*delta2;      % Eq 30    
end

vars = [mu2(:) mu3(:) sigma2(:) sigma3(:)];
if any(isnan(vars(:)))
    bad_traj = true;
end

m = mu1hat(2:end,:);
end
