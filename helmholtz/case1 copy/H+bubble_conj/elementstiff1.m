function [K] = elementstiff1(X, Y, m, n, k0)
%stiffness element on the local patch to find basis on each coarse mesh 
K = zeros(4, 4);
vs = cellVertices(X, Y, m, n);
xlow = vs(1, 1);
xhigh = vs(3, 1);
ylow = vs(1, 2);
yhigh = vs(3, 2);

x = (xlow+xhigh)/2;
y = (ylow+yhigh)/2;
length = -xlow+xhigh;
beta = [betafun(x,ylow)+betafun(xlow,y),betafun(x,ylow)+betafun(xhigh,y),betafun(x,yhigh)+betafun(xhigh,y),betafun(x,yhigh)+betafun(xlow,y)];
beta1 = [betafun(x,ylow),betafun(xhigh,y),betafun(x,yhigh)];
for i = 1:4 
    for j = 1:4
        if i==j
            K(i,j)=2/3*afun(x,y)-k0^2*vfun(x,y)^2*1/9*length^2-1i*k0*beta(i)*1/3*length;
        elseif i==j+2 || i==j-2
                K(i,j)=-1/3*afun(x,y)-k0^2*vfun(x,y)^2*1/36*length^2;
        elseif i==j+1 || i==j-1
            index = min(i,j);
            K(i,j)=-1/6*afun(x,y)-k0^2*vfun(x,y)^2*1/18*length^2-1i*k0*beta1(index)*1/6*length;
        else
            K(i,j)=-1/6*afun(x,y)-k0^2*vfun(x,y)^2*1/18*length^2-1i*k0*betafun(xlow,y)*1/6*length;
        end
    end
end
end





