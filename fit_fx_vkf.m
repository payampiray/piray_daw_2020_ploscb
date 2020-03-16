function [output] = fit_fx_vkf(params)
lux      = @(x)(1/(1+exp(-x)));
lambda   = lux(params(1));
v0       = 10*lux(params(2));
omega    = exp(params(3));
beta     = exp(params(4));

output = struct('lambda',lambda,'v0',v0,'omega',omega,'beta',beta);

end