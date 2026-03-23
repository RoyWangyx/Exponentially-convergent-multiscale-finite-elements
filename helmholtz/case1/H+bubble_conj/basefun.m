function [A] = basefun(X, Y, m, n,N_f,k0)
% assemble the local stiffness matrix for each coarse element based on 
% the stiffness elements computed in elementstiff1

x = linspace(X(m), X(m+1), N_f+1);
y = linspace(Y(n), Y(n+1), N_f+1);
nNodes = (N_f+1)*(N_f+1);
%sparse assembling
I=zeros(16*N_f^2,1);J=I;K=I;

for i = 1:N_f
    for j = 1:N_f  
        [k] = elementstiff1(x, y, i, j,k0);
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

A=sparse(I,J,K,nNodes,nNodes);


end






