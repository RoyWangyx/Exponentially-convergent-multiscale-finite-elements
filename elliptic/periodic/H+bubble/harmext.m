function [L1,L2,N] = harmext(X, Y, m, n,N_f,i)
%compute the local harmonic extension corresponding to each coarse edge
%i=1, corresponds to the horizontal edges, i=2, corresponds to the vertical
%ones. m, n are indices

%use dirichlet solver for adajacent patches, obtaining L1.L2 as a linear
%combination of the corresponding patches, N is matrix of the inner product

b=[1:N_f+1,N_f+2:N_f+1:(N_f+1)*(N_f+1),2*(N_f+1):N_f+1:(N_f+1)*(N_f+1),N_f*N_f+N_f+2:(N_f+1)*(N_f+1)-1];
K1 = basefun(X, Y, m, n,N_f);
f1= zeros(4*N_f,N_f-1);
f2= zeros(4*N_f,N_f-1);
if i==1
    K2 = basefun(X, Y, m, n+1,N_f);
    f1(3*N_f+2:4*N_f,:)= speye(N_f-1);
    f2(2:N_f,:)= speye(N_f-1);
else
    K2 = basefun(X, Y, m+1, n,N_f);
    f1(2*N_f+2:3*N_f,:)= speye(N_f-1);
    f2(N_f+2:2*N_f,:)= speye(N_f-1);
end
M1=K1;
M2=K2;
F1=-K1(:,b)*f1;
K1(b,:)=0; K1(:,b)=0; F1(b,:)=f1; K1(b,b)=speye(length(b),length(b));
L1=K1\F1;
F2=-K2(:,b)*f2;
K2(b,:)=0; K2(:,b)=0; F2(b,:)=f2; K2(b,b)=speye(length(b),length(b));
L2=K2\F2;

N1=L1'*M1*L1;
N2=L2'*M2*L2;
N=N1+N2;
   


end
