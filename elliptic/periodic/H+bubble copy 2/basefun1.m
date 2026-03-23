function [A,N_x,N_y] = basefun1(X, Y, m, n,N_f,t)
% solve a combination of fine basis functions, assembling the oversampled
% stiffness matrix
%t=1, corresponds to the horizontal edges, t=2, corresponds to the vertical
%ones. m, n are indices

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
N_c=length(X)-1;
if t==1
    if m==1
        x = linspace(X(m), X(m+2), 2*N_f+1);
    elseif m==N_c
        x = linspace(X(m-1), X(m+1), 2*N_f+1);
    else
        x = linspace(X(m-1), X(m+2), 3*N_f+1);
    end
    y = linspace(Y(n), Y(n+2), 2*N_f+1);
else
    if n==1
        y = linspace(Y(n), Y(n+2), 2*N_f+1);
    elseif n==N_c
        y = linspace(Y(n-1), Y(n+1), 2*N_f+1);
    else
        y = linspace(Y(n-1), Y(n+2), 3*N_f+1);
    end
    x = linspace(X(m), X(m+2), 2*N_f+1);
end
    
    


nNodes = length(x)*length(y);
N_x=length(x)-1;
N_y=length(y)-1;
%sparse assembling
I=zeros(16*N_x*N_y,1);J=I;K=I;

for i = 1:N_x
    for j = 1:N_y  
        [k] = elementstiff1(x, y, i, j);
%         idx = loc2glo(N, m, n, i)
        for p = 1:4
            global_p = loc2glo(N_x, i, j, p);
            for q = 1:4
                index=16*N_y*(i-1)+16*(j-1)+4*(p-1)+q;
                global_q = loc2glo(N_x, i, j, q);
                I(index)=global_p;
                J(index)=global_q;
                K(index)=k(p, q);
            end
        end
    end
end

%assemble boundary

A=sparse(I,J,K,nNodes,nNodes);


end






