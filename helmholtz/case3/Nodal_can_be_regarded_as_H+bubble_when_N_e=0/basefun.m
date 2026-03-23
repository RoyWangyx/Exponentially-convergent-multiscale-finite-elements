function [B,value] = basefun(X, Y, m, n,N_f)
% solve a combination of fine basis functions

% the function serves as to solve local boundary value problem on the 
% coarse mesh, based on the stiffness matrix in "elementstiff1"

% for the sake of simplicity, we temporarily use nodes of the fine mesh
% as the same nodes that we use when we invoke numerical quadrature to
% assemble the global stiffness matrix, the rationale being that we
% can either increase the scale of fine mesh or increase the quadrature 
% accuracy if we want to do higher order.... (because locally on fine
% meshes, the basis function is simply a linear function

% setting parameters H is the size of the coarse mesh, N is the number of
%nodes on each coarse rectangular, i.e. fine mesh size is H/N


%N needs to be inputted !!!!

x = linspace(X(m), X(m+1), N_f+1);
y = linspace(Y(n), Y(n+1), N_f+1);
nNodes = (N_f+1)*(N_f+1);
%sparse assembling
I=zeros(16*N_f^2,1);J=I;K=I;

for i = 1:N_f
    for j = 1:N_f  
        [k] = elementstiff1(x, y, i, j);
%         idx = loc2glo(N, m, n, i)
        for p = 1:4
            global_p = loc2glo(N_f, i, j, p);
            for q = 1:4
                index=16*N_f*(i-1)+16*(j-1)+4*(p-1)+q;
                global_q = loc2glo(N_f, i, j, q);
                I(index)=global_p;
                J(index)=global_q;
                K(index)=k(p, q);
            end
        end
    end
end

%assemble boundary
b=[1:N_f+1,N_f+2:N_f+1:(N_f+1)*(N_f+1),2*(N_f+1):N_f+1:(N_f+1)*(N_f+1),N_f*N_f+N_f+2:(N_f+1)*(N_f+1)-1];
f=[linspace(1, 0, N_f+1), linspace(1-1/N_f, 0, N_f),zeros(1,N_f),zeros(1,N_f-1);...
    linspace(0, 1, N_f+1),zeros(1,N_f),linspace(1-1/N_f, 0, N_f),zeros(1,N_f-1);...
    zeros(1,N_f+1),zeros(1,N_f),linspace(1/N_f, 1, N_f),linspace(1/N_f, 1-1/N_f, N_f-1);...
    zeros(1,N_f+1),linspace(1/N_f, 1, N_f),zeros(1,N_f),linspace(1-1/N_f, 1/N_f, N_f-1)];
A=sparse(I,J,K,nNodes,nNodes);
F=-A(:,b)*f';
B=A;
A(b,:)=0; A(:,b)=0; F(b,:)=f'; A(b,b)=speye(length(b),length(b));
value=A\F;

end






