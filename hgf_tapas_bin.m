function [tapas_m, bad_traj, tapas_mu2, tapas_mu3, tapas_sigma2, e] = hgf_tapas_bin(o,nu,kappa,omega)
lognu = log(nu);
p = [NaN, 0, 1.0000, NaN, 0.1000, 1.0000, NaN, 0, 0, 1.0000, kappa, NaN, omega, lognu];

r.y = ones(size(o));
r.u = o;
r.ign = [];
r.irr = zeros(0,1);
r.c_prc = tapas_hgf_binary_config;
r.c_obs = tapas_unitsq_sgm_config;
r.c_opt = tapas_quasinewton_optim_config;

bad_traj = 0;
tapas_m = nan;
tapas_mu2 = nan;
tapas_mu3 = nan;
tapas_sigma2 = nan;
e = nan;
try
    [traj] = tapas_hgf_binary(r, p);
catch
    bad_traj = 1;
    return;
end

tapas_m = traj.muhat(:,1);
tapas_mu2 = [0; traj.mu(:,2)];
tapas_mu3 = [1; traj.mu(:,3)];
tapas_sigma2 = [0.1; traj.sa(:,2)];

%-----------
% compare the output with my implementation
[m, ~, mu2, mu3, sigma2] =hgf_bin(o,nu,kappa,omega);

e1 = max(abs(m - tapas_m));
e2 = max(abs(mu2 - tapas_mu2));
e3 = max(abs(mu3 - tapas_mu3));
e4 = max(abs(sigma2 - tapas_sigma2));
e = [e1 e2 e3 e4];

end
