function f=ffun(x,y)
if (x-0.5)^2+(y-0.5)^2<1/400
    f = exp(-1/(1-400*((x-0.5)^2+(y-0.5)^2)));
else
    f= 0;
end
end