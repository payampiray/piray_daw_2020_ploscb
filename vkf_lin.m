function [val, lr, vol] = vkf_lin(y,lambda,v0,sigma2)

[nt,nq] = size(y);

m       = zeros(1,nq);
w       = sigma2*ones(1,nq);
v       = v0*ones(1,nq);

val     = nan(nt,nq);
lr      = nan(nt,nq);
vol     = nan(nt,nq);

for t  = 1:nt              
    
    wpre        = w;            
    
    delta       = y(t,:) - m;
    k           = (w+v)./(w+v + sigma2);    
    m           = m + k.*delta;    
    w           = (1-k).*(w+v);
    
    v           = v +lambda.*(k.^2.*delta.^2 + k.*wpre - k.*v);
    
    val(t,:)    = m;
    vol(t,:)    = v;        
    lr(t,:)     = k;    
end

end
