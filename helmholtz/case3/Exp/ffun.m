function f=ffun(x,y)
if (x-0.125)^2+(y-0.5)^2<1/400
    f = 10000*exp(-1/(1-400*((x-0.125)^2+(y-0.5)^2)));
else
    f= 0;
end
end