function [m_signal,v_signal,cc,pf] = rbpf_hgf_lin(y,nu,kappa,omega,alpha,nparticles)
z_rng = [1 1];
% x_rng = [0 0];

state_model = @(particles)pf_state_transition(particles, nu);
measurement_model = @pf_measurement_likelihood;

pf = particleFilter(state_model,measurement_model);
initialize(pf,nparticles,z_rng);

pf.StateEstimationMethod = 'mean';
pf.ResamplingMethod = 'systematic';

m = zeros(1,nparticles);
w = ones(1,nparticles);

N = length(y);
xEst = zeros(N,2);
cc(1,:) = zeros(1,1);
for t=1:size(y,1)
    xEst(t,1) = predict(pf);
    xEst(t,2) = pf.Weights*m';    
    correct(pf,y(t),m,w,kappa,omega,alpha);
    [m,w]=kalman(pf.Particles,y(t),m,w,kappa,omega,alpha);
    [~,ccov]=getStateEstimate(pf);
    cc(t+1,:) = ccov;
end
xEst(t+1,1) = predict(pf);
xEst(t+1,2) = pf.Weights*m';

v_signal = xEst(:,1);
m_signal = xEst(:,2);

end


%------------------------------
function particles = pf_state_transition(particles, nu)
% number of states x number of particles

% the mean of particles does not change for normal state-space model
z = particles(1,:);
z = z + sqrt(nu).*randn(size(z));
particles = z;

end

function likelihood = pf_measurement_likelihood(particles, measurement, m,w,kappa,omega,alpha)
% the mean of particles does not change for normal state-space model

z = particles;
sigma = exp(kappa*z+omega);

% x = particles(2,:);
likelihood = normpdf(measurement,m,sqrt(w+sigma+alpha));

end

function [m,w]=kalman(particles,yt,m,w,kappa,omega,alpha)
z = particles;
sigma = exp(kappa*z+omega);

k = (w+sigma)./(w+sigma+alpha);
m = m + k.*(yt-m);
w = (1-k).*(w+sigma);
end
