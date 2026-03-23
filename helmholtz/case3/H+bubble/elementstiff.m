function [value,K, f] = elementstiff(X, Y, m, n,N_f,N_e,k0)
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
if m == 1 && 1<n && n<N_c
    b=[1:N_f+1,2*(N_f+1):N_f+1:(N_f+1)*(N_f+1),N_f*N_f+N_f+1:(N_f+1)*(N_f+1)-1];
    f=[linspace(1, 0, N_f+1), zeros(1,N_f),zeros(1,N_f);...
    linspace(0, 1, N_f+1),linspace(1-1/N_f, 0, N_f),zeros(1,N_f);...
    zeros(1,N_f+1),linspace(1/N_f, 1, N_f),linspace(0, 1-1/N_f, N_f);...
    zeros(1,N_f+1),zeros(1,N_f),linspace(1, 1/N_f, N_f)];
elseif m == N_c && 1<n && n<N_c
    b=[1:N_f+1,N_f+2:N_f+1:(N_f+1)*(N_f+1),N_f*N_f+N_f+2:(N_f+1)*(N_f+1)];
    f=[linspace(1, 0, N_f+1), linspace(1-1/N_f, 0, N_f),zeros(1,N_f);...
    linspace(0, 1, N_f+1),zeros(1,N_f),zeros(1,N_f);...
    zeros(1,N_f+1),zeros(1,N_f),linspace(1/N_f, 1, N_f);...
    zeros(1,N_f+1),linspace(1/N_f, 1, N_f),linspace(1-1/N_f, 0, N_f)];
elseif n == 1 && 1<m && m<N_c
    b=[1:N_f+1:(N_f+1)*(N_f+1),(N_f+1):N_f+1:(N_f+1)*(N_f+1),N_f*N_f+N_f+2:(N_f+1)*(N_f+1)-1];
    f=[linspace(1, 0, N_f+1),zeros(1,N_f+1),zeros(1,N_f-1);...
    zeros(1,N_f+1),linspace(1, 0, N_f+1),zeros(1,N_f-1);...
    zeros(1,N_f+1),linspace(0, 1, N_f+1),linspace(1/N_f, 1-1/N_f, N_f-1);...
    linspace(0, 1, N_f+1),zeros(1,N_f+1),linspace(1-1/N_f, 1/N_f, N_f-1)];
elseif n == N_c && 1<m && m<N_c
    b=[1:N_f+1,N_f+2:N_f+1:(N_f+1)*(N_f+1),2*(N_f+1):N_f+1:(N_f+1)*(N_f+1)];
    f=[linspace(1, 0, N_f+1), linspace(1-1/N_f, 0, N_f),zeros(1,N_f);...
    linspace(0, 1, N_f+1),zeros(1,N_f),linspace(1-1/N_f, 0, N_f);...
    zeros(1,N_f+1),zeros(1,N_f),linspace(1/N_f, 1, N_f);...
    zeros(1,N_f+1),linspace(1/N_f, 1, N_f),zeros(1,N_f)];
elseif m == 1 && n ==1
    b=[(N_f+1):N_f+1:(N_f+1)*(N_f+1),N_f*N_f+N_f+1:(N_f+1)*(N_f+1)-1];
    f=[zeros(1,N_f+1),zeros(1,N_f);...
    linspace(1, 0, N_f+1),zeros(1,N_f);...
    linspace(0, 1, N_f+1),linspace(0, 1-1/N_f, N_f);...
    zeros(1,N_f+1),linspace(1, 1/N_f, N_f)];
elseif m == 1 && n ==N_c
    b=[1:N_f+1,2*(N_f+1):N_f+1:(N_f+1)*(N_f+1)];
    f=[linspace(1, 0, N_f+1),zeros(1,N_f);...
    linspace(0, 1, N_f+1),linspace(1-1/N_f, 0, N_f);...
    zeros(1,N_f+1),linspace(1/N_f, 1, N_f);...
    zeros(1,N_f+1),zeros(1,N_f)];
elseif m == N_c && n ==1
    b=[1:N_f+1:(N_f+1)*(N_f+1),N_f*N_f+N_f+2:(N_f+1)*(N_f+1)];
    f=[linspace(1, 0, N_f+1),zeros(1,N_f);...
    zeros(1,N_f+1),zeros(1,N_f);...
    zeros(1,N_f+1),linspace(1/N_f, 1, N_f);...
    linspace(0, 1, N_f+1),linspace(1-1/N_f, 0, N_f)];
elseif m == N_c && n ==N_c
    b=[1:N_f+1,N_f+2:N_f+1:(N_f+1)*(N_f+1)];
    f=[linspace(1, 0, N_f+1), linspace(1-1/N_f, 0, N_f);...
    linspace(0, 1, N_f+1),zeros(1,N_f);...
    zeros(1,N_f+1),zeros(1,N_f);...
    zeros(1,N_f+1),linspace(1/N_f, 1, N_f)];
else
    b=[1:N_f+1,N_f+2:N_f+1:(N_f+1)*(N_f+1),2*(N_f+1):N_f+1:(N_f+1)*(N_f+1),N_f*N_f+N_f+2:(N_f+1)*(N_f+1)-1];
    f=[linspace(1, 0, N_f+1), linspace(1-1/N_f, 0, N_f),zeros(1,N_f),zeros(1,N_f-1);...
    linspace(0, 1, N_f+1),zeros(1,N_f),linspace(1-1/N_f, 0, N_f),zeros(1,N_f-1);...
    zeros(1,N_f+1),zeros(1,N_f),linspace(1/N_f, 1, N_f),linspace(1/N_f, 1-1/N_f, N_f-1);...
    zeros(1,N_f+1),linspace(1/N_f, 1, N_f),zeros(1,N_f),linspace(1-1/N_f, 1/N_f, N_f-1)];
end
A=basefun(X, Y, m, n,N_f,k0);
B=A;
F=-A(:,b)*f';
A(b,:)=0; A(:,b)=0; F(b,:)=f'; A(b,b)=speye(length(b),length(b));
value=A\F;
%edge basis
count=4;
if n>1
    [L1,L2,N]=harmext(X, Y, m, n-1,N_f,1,k0);
    [R,P] = restrict(X, Y, m, n-1,N_f,1,k0);
    [V,D]=eigs(R'*N*R,P,N_e);
    value(:,count+1:count+N_e)=L2*R*V;
else
    value(:,count+1:count+N_e)=0;
end
count=count+N_e;
if n<N_c
    [L1,L2,N]=harmext(X, Y, m, n,N_f,1,k0);
    [R,P] = restrict(X, Y, m, n,N_f,1,k0);
    [V,D]=eigs(R'*N*R,P,N_e);
    value(:,count+1:count+N_e)=L1*R*V;
else
    value(:,count+1:count+N_e)=0;
end
count=count+N_e;
if m>1
    [L1,L2,N]=harmext(X, Y, m-1, n,N_f,2,k0);
    [R,P] = restrict(X, Y, m-1, n,N_f,2,k0);
    [V,D]=eigs(R'*N*R,P,N_e);
    value(:,count+1:count+N_e)=L2*R*V;
else
    value(:,count+1:count+N_e)=0;
end
count=count+N_e;
if m<N_c
    [L1,L2,N]=harmext(X, Y, m, n,N_f,2,k0);
    [R,P] = restrict(X, Y, m, n,N_f,2,k0);
    [V,D]=eigs(R'*N*R,P,N_e);
    value(:,count+1:count+N_e)=L1*R*V;
else
    value(:,count+1:count+N_e)=0;
end
count=count+N_e;

bubble_tmp=reshape(bubble(X,Y,m,n,N_f,k0),(N_f+1)^2,1);

K=value.'*B*conj(value);
%k and f separate using local stiffness matrix%%%%%%
F = zeros(N_f^2+4*N_f,count);
for u = 1:N_f
    xmean=(x(u+1)+x(u))/2;
    for v = 1:N_f
        ymean=(y(v+1)+y(v))/2;
        for x1=0:1
            for y1=0:1
                F(v*N_f+u-N_f,:) = F(v*N_f+u-N_f,:)+ffun(xmean,ymean)*conj(value((v+y1-1)*(N_f+1)+u+x1,:))*(x(2)-x(1))^2/4;
            end
        end
    end
end
for u = 1:N_f
    xmean=(x(u+1)+x(u))/2;
    ymean=(y(u+1)+y(u))/2;
    for x1=0:1
        F(N_f*N_f+u,:) = F(N_f*N_f+u,:)+gfun(xmean,yhigh,k0)*conj(value(N_f*(N_f+1)+u+x1,:))*(x(2)-x(1))/2;
        F(N_f*N_f+u+N_f,:) = F(N_f*N_f+u+N_f,:)+gfun(xmean,ylow,k0)*conj(value(u+x1,:))*(x(2)-x(1))/2;
        F(N_f*N_f+u+2*N_f,:) = F(N_f*N_f+u+2*N_f,:)+gfun(xhigh,ymean,k0)*conj(value((u+x1-1)*(N_f+1)+N_f+1,:))*(x(2)-x(1))/2;
        F(N_f*N_f+u+3*N_f,:) = F(N_f*N_f+u+3*N_f,:)+gfun(xlow,ymean,k0)*conj(value((u+x1-1)*(N_f+1)+1,:))*(x(2)-x(1))/2;
    end
end
f=sum(F);
f=f-bubble_tmp.'*B*conj(value);
end


