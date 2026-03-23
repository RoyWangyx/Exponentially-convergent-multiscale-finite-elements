function [K, f] = elementstiff(X, Y, m, n,N_f)
%compute the global stiffness matrix

vs = cellVertices(X, Y, m, n);
xlow = vs(1, 1);
xhigh = vs(3, 1);

x = linspace(xlow, xhigh, N_f+1);


[B,value] = basefun(X, Y, m, n,N_f);
K=value'*B*value;
%k and f separate using local stiffness matrix%%%%%%
f = zeros(4, 1);
for i = 1:4
    F = zeros(N_f, N_f);
    for u = 1:N_f
        for v = 1:N_f
            for x1=0:1
                for y1=0:1
                    F(u, v) = F(u,v)-value((v+y1-1)*(N_f+1)+u+x1,i)*(x(2)-x(1))^2/4;
                end
            end
        end
    end
    f(i)=sum(sum(F));
end
end


