function [L1,L2,N] = harmext(X, Y, m, n,N_f,t,k0)
%compute the local harmonic extension corresponding to each coarse edge
%i=1, corresponds to the horizontal edges, i=2, corresponds to the vertical
%ones. m, n are indices

%use dirichlet solver for adjacent patches, obtaining L1.L2 as a linear
%combination of the corresponding patches, N is matrix of the inner product


N_c=length(X)-1;
b1 = bc(N_c, N_c,N_f, N_f, m, n);
D1 = size(b1,2);
K1 = basefun(X, Y, m, n,N_f,k0);
f1= zeros(D1,N_f-1);
if t==1
    K2 = basefun(X, Y, m, n+1,N_f,k0);
    b2 = bc(N_c, N_c,N_f, N_f,m, n+1);
    D2 = size(b2,2);
    f2= zeros(D2,N_f-1);
    [tf1,loc1] = ismember([N_f*N_f+N_f+2:(N_f+1)*(N_f+1)-1],b1);
    [tf2,loc2] = ismember([2:N_f],b2);
    f1(loc1,:)= speye(N_f-1);
    f2(loc2,:)= speye(N_f-1);
else
    K2 = basefun(X, Y, m+1, n,N_f,k0);
    b2 = bc(N_c, N_c,N_f, N_f,m+1, n);
    D2 = size(b2,2);
    f2= zeros(D2,N_f-1);
    [tf2,loc2] = ismember([N_f+2:N_f+1:(N_f+1)*N_f],b2);
    [tf1,loc1] = ismember([2*(N_f+1):N_f+1:(N_f+1)*(N_f)],b1);
    f1(loc1,:)= speye(N_f-1);
    f2(loc2,:)= speye(N_f-1);
end
M1=K1;
M2=K2;
F1=-K1(:,b1)*f1;
K1(b1,:)=0; K1(:,b1)=0; F1(b1,:)=f1; K1(b1,b1)=speye(length(b1),length(b1));
L1=K1\F1;
F2=-K2(:,b2)*f2;
K2(b2,:)=0; K2(:,b2)=0; F2(b2,:)=f2; K2(b2,b2)=speye(length(b2),length(b2));
L2=K2\F2;

N1=L1'*M1*L1;
N2=L2'*M2*L2;
N=N1+N2;
   


end
