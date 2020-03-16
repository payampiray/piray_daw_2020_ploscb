function [m_signal,v_signal,pf] = pf_vkf_bin(y,lambda,v0,nparticles)
z_rng = [v0 v0].^-1;
x_rng = [0 0];

state_model = @(particles)pf_state_transition(particles, lambda);
measurement_model = @(particles, measurement)pf_measurement_likelihood(particles, measurement);

pf = particleFilter(state_model,measurement_model);
initialize(pf,nparticles,[z_rng;x_rng]);

pf.StateEstimationMethod = 'mean';
pf.ResamplingMethod = 'systematic';

N = length(y);
xEst = zeros(N,2);
for t=1:size(y,1)
    xEst(t,:) = predict(pf);
    correct(pf,y(t));    
end
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

x = particles(2,:);
x = x + sqrt(z.^-1).*randn(size(x));

particles = [z; x];
end

function likelihood = pf_measurement_likelihood(particles, measurement)
% the mean of particles does not change for normal state-space model

x = particles(2,:);
p = 1./(1+exp(-x));
% p = normcdf(x);
o = measurement;
likelihood = p.*o + (1-p).*(1-o);
end