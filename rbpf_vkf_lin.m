function [m_signal,v_signal,cc,pf] = rbpf_vkf_lin(y,lambda,v0,omega,nparticles)
z_rng = [v0 v0].^-1;
% x_rng = [0 0];

state_model = @(particles)pf_state_transition(particles, lambda);
measurement_model = @(particles, measurement)pf_measurement_likelihood(particles, measurement, omega);
measurement_model = @pf_measurement_likelihood;

pf = particleFilter(state_model,measurement_model);
initialize(pf,nparticles,z_rng);

pf.StateEstimationMethod = 'mean';
pf.ResamplingMethod = 'systematic';

N = length(y);
xEst = zeros(N,2);
cc(1,:) = zeros(1,4);

m = zeros(1,nparticles);
w = ones(1,nparticles);

for t=1:size(y,1)    
    xEst(t,1) = predict(pf);
    xEst(t,2) = pf.Weights*m';
    correct(pf,y(t),m,w,omega);
    [m,w]=kalman(pf.Particles,y(t),m,w,omega);    
    [s,ccov]=getStateEstimate(pf);
    cc(t+1,:) = ccov;    
end
xEst(t+1,1) = predict(pf);
xEst(t+1,2) = pf.Weights*m';

v_signal = xEst(:,1).^-1;
m_signal = xEst(:,2);

end


%------------------------------
function particles = pf_state_transition(particles, lambda)
% number of states x number of particles

% the mean of particles does not change for normal state-space model
z = particles(1,:);
eta = 1-lambda;
nu = .5/(1-eta);
epsil = betarnd(eta*nu,(1-eta)*nu, size(z));
e = (eta.^-1)*epsil;
z = z.*e;
particles = z;

end

function likelihood = pf_measurement_likelihood(particles, measurement, m,w,omega)
% the mean of particles does not change for normal state-space model

z = particles;
sigma = z.^-1;

likelihood = normpdf(measurement,m,sqrt(w+sigma+omega));

end

function [m,w]=kalman(particles,yt,m,w,omega)
z = particles;
sigma = z.^-1;

k = (w+sigma)./(w+sigma+omega);
m = m + k.*(yt-m);
w = (1-k).*(w+sigma);
end