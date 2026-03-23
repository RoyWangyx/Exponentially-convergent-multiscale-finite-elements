function [A,N_x,N_y,x,y] = basefun1(X, Y, m, n,N_f,t,k0)
% assemble the oversampling local stiffness matrix for each coarse element based on 
% the stiffness elements computed in elementstiff1
%t=1, corresponds to the horizontal edges, t=2, corresponds to the vertical

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
        [k] = elementstiff1(x, y, i, j,k0);
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






