function [f] =gfun(x,y,k0)
a = [0.2, 1.6, 1.8, 0.4, 0];
if x==0
    j = 4;
elseif x==1
    j = 2;
elseif y ==0
    j = 1;
elseif y==1
    j = 3;
else
    j = 5;
end
f = -1i*k0*a(j)*exp(-0.6i*k0*x-0.8i*k0*y);
end