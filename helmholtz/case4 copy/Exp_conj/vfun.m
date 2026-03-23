function [f] =vfun(t,s)
global gov;
M=128;
t1=floor(M*t)+1;
t2=floor(M*t)+2;
s1=floor(M*s)+1;
s2=floor(M*s)+2;
f=(s1-M*s)*(t1-M*t)*gov(t1,s1)+(s1-M*s)*(M*t+1-t1)*gov(t2,s1)+...
    (M*s+1-s1)*(t1-M*t)*gov(t1,s2)+(M*s+1-s1)*(M*t+1-t1)*gov(t2,s2);
end