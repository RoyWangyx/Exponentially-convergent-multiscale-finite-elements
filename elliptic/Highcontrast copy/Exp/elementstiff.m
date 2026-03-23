function [value,K, f] = elementstiff(X, Y, m, n,N_f,N_e)
%compute the global stiffness matrix
N_c=length(X)-1;
vs = cellVertices(X, Y, m, n);
xlow = vs(1, 1);
xhigh = vs(3, 1);
ylow = vs(1, 2);
yhigh = vs(3, 2);
x = linspace(xlow, xhigh, N_f+1);
y = linspace(ylow, yhigh, N_f+1);
%nodal basis
b=[1:N_f+1,N_f+2:N_f+1:(N_f+1)*(N_f+1),2*(N_f+1):N_f+1:(N_f+1)*(N_f+1),N_f*N_f+N_f+2:(N_f+1)*(N_f+1)-1];
f=[linspace(1, 0, N_f+1), linspace(1-1/N_f, 0, N_f),zeros(1,N_f),zeros(1,N_f-1);...
    linspace(0, 1, N_f+1),zeros(1,N_f),linspace(1-1/N_f, 0, N_f),zeros(1,N_f-1);...
    zeros(1,N_f+1),zeros(1,N_f),linspace(1/N_f, 1, N_f),linspace(1/N_f, 1-1/N_f, N_f-1);...
    zeros(1,N_f+1),linspace(1/N_f, 1, N_f),zeros(1,N_f),linspace(1-1/N_f, 1/N_f, N_f-1)];
A=basefun(X, Y, m, n,N_f);
B=A;
F=-A(:,b)*f';
A(b,:)=0; A(:,b)=0; F(b,:)=f'; A(b,b)=speye(length(b),length(b));
value=A\F;
%edge basis
count=4;
if n>1
    [L1,L2,N]=harmext(X, Y, m, n-1,N_f,1);
    [R,P,bub] = restrict(X, Y, m, n-1,N_f,1);
    [V,D]=eigs(R'*N*R,P,N_e);
    value(:,count+1:count+N_e)=L2*R*V;
    value(:,count+N_e+1)=L2*bub;
else
    value(:,count+1:count+N_e+1)=0;
end
count=count+N_e+1;
if n<N_c
    [L1,L2,N]=harmext(X, Y, m, n,N_f,1);
    [R,P,bub] = restrict(X, Y, m, n,N_f,1);
    [V,D]=eigs(R'*N*R,P,N_e);
    value(:,count+1:count+N_e)=L1*R*V;
    value(:,count+N_e+1)=L1*bub;
else
    value(:,count+1:count+N_e+1)=0;
end
count=count+N_e+1;
if m>1
    [L1,L2,N]=harmext(X, Y, m-1, n,N_f,2);
    [R,P,bub] = restrict(X, Y, m-1, n,N_f,2);
    [V,D]=eigs(R'*N*R,P,N_e);
    value(:,count+1:count+N_e)=L2*R*V;
    value(:,count+N_e+1)=L2*bub;
else
    value(:,count+1:count+N_e+1)=0;
end
count=count+N_e+1;
if m<N_c
    [L1,L2,N]=harmext(X, Y, m, n,N_f,2);
    [R,P,bub] = restrict(X, Y, m, n,N_f,2);
    [V,D]=eigs(R'*N*R,P,N_e);
    value(:,count+1:count+N_e)=L1*R*V;
    value(:,count+N_e+1)=L1*bub;
else
    value(:,count+1:count+N_e+1)=0;
end
count=count+N_e+1;



K=value'*B*value;
%k and f separate using local stiffness matrix%%%%%%
F = zeros(N_f^2,count);
for u = 1:N_f
    for v = 1:N_f
        for x1=0:1
            for y1=0:1
                xf=(x(u+1)+x(u))/2;
                yf=(y(v+1)+y(v))/2;
                F(v*N_f+u-N_f,:) = F(v*N_f+u-N_f,:)+ffun(xf,yf)*value((v+y1-1)*(N_f+1)+u+x1,:)*(x(2)-x(1))^2/4;
            end
        end
    end
end
f=sum(F);
end


