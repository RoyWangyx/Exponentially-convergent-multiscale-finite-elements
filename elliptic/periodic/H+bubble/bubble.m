function [re] = bubble(X, Y, m, n,N_f)
%solve the bubble part
x = linspace(X(m), X(m+1), N_f+1);
y = linspace(Y(n), Y(n+1), N_f+1);
nNodes = (N_f+1)*(N_f+1);
%sparse assembling
I=zeros(16*N_f^2,1);J=I;K=I;
F=zeros(nNodes, 1);

for i = 1:N_f
    for j = 1:N_f  
        [k,f] = elementstiff2(x, y, i, j);
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
            F(global_p) = F(global_p) + f(p);
        end
    end
end



  
%assemble boundary
b=[1:N_f+1,N_f+2:N_f+1:(N_f+1)*(N_f+1),2*(N_f+1):N_f+1:(N_f+1)*(N_f+1),N_f*N_f+N_f+2:(N_f+1)*(N_f+1)-1,];

A=sparse(I,J,K,nNodes,nNodes);

A(b,:)=0; A(:,b)=0; F(b)=0; A(b,b)=speye(length(b),length(b));
u = A\F;
re=reshape(u,N_f+1,N_f+1);

end