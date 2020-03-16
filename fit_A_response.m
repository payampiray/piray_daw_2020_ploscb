function [loglik, beta] = fit_A_response(X,Y,params)
beta     = exp(params(1));
bv       = params(2);
be       = params(3);
bi       = params(4);
bb       = [bv bv -bv -bv] + [be -be be -be] + [bi -bi -bi bi];

Y        = logical(Y);
z        = bsxfun(@plus,X*beta , bb);
f        = (1./(1+exp(-z)));
    
p        = f.*Y + (1-f).*(1-Y);
loglik   = sum(sum(log(p+eps)));
end