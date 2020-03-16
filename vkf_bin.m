function [dv, lr, vol, um] =vkf_bin(y,lambda,v0,omega)

[nt,nq] = size(y);

m       = 0*ones(1,nq);
w       = omega*ones(1,nq);
v       = v0*ones(1,nq);

dv      = nan(nt,nq);
lr      = nan(nt,nq);
vol     = nan(nt,nq);
um      = nan(nt,nq);


ux = @(x)1./(1+exp(-x));
for t  = 1:nt      
    dv(t,:)     = m(1,:);
    lr(t,:)     = sqrt(w+v);
    vol(t,:)    = v;
    
    wpre        = w;            
    
    delta       = sqrt(w+v).*(y(t,:) - ux(m));
    m           = m + delta;

    k           = (w+v)./(w+v + omega);
    w           = (1-k).*(w+v);
    
    v           = v +lambda.*(delta.^2 + k.*wpre - k.*v);
    
    um(t,:)     = ux(m);
end

end
