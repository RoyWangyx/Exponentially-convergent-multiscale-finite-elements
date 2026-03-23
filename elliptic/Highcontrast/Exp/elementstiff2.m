function [K,f] = elementstiff2(X, Y, m, n)
%the stiffness matrix required locally to compute local nodal basis 
%for the bilinear boundary value with 1 at node "node"

%X Y are the fine mesh, we shall use local bilinear basis on the 
%fine mesh as basis to compute the harmonic problem for bilinear
%boundary value problem on the coarse mesh

%H is the size of the coarse mesh 


%H, nSamples and coefficients a needs to be inputted  !!!!!!!!!

%nSamples = 3; % sample points in trapz in each direction
K = zeros(4, 4);
vs = cellVertices(X, Y, m, n);
xlow = vs(1, 1);
xhigh = vs(3, 1);
ylow = vs(1, 2);
yhigh = vs(3, 2);
f=zeros(4,1);
x = (xlow+xhigh)/2;
y = (ylow+yhigh)/2;

for i = 1:4 
    for j = 1:4
        if i==j
            K(i, j)=2/3*afun(x,y);
        elseif i==j+2 || i==j-2
                K(i,j)=-1/3*afun(x,y);
        else 
            K(i,j)=-1/6*afun(x,y);
        end
    end
end

for i=1:4
    f(i) = ffun(x,y)*(xhigh-xlow)^2/4;
end

end