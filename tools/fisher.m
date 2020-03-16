function z  = fisher(x)
x(x==-1) = x(x==-1)+eps;
x(x==+1) = x(x==+1)-eps;
z = .5*log( (1+x)./(1-x) );

end