clc; clear;
tic
% setting parameters
N_c= 16;

N_f=16;
nSamples = N_f+1;
fprintf('Coarse Grid %d X %d \n', N_c, N_c);
fprintf('Fine Grid %d X %d \n', N_c*N_f, N_c*N_f);
usebilinear = true;
% domain [0, 1]X[0, 1]
x = linspace(0, 1, N_c+1);
y = x;
nNodes = (N_c+1)*(N_c+1);

%sparse assembling
I=zeros(16*N_c^2,1);J=I;K=I;
F = zeros(nNodes, 1);


for i = 1:N_c
    for j = 1:N_c  
        [k, f] = elementstiff(x, y, i, j,N_f);
%         idx = loc2glo(N, m, n, i)
        for p = 1:4
            global_p = loc2glo(N_c, i, j, p);
            for q = 1:4
                index=16*N_c*(i-1)+16*(j-1)+4*(p-1)+q;
                global_q = loc2glo(N_c, i, j, q);
                I(index)=global_p;
                J(index)=global_q;
                K(index)=k(p, q);
            end
            F(global_p) = F(global_p) + f(p);
        end
    end
end

b=[1:N_c+1,N_c+2:N_c+1:(N_c+1)*(N_c+1),2*(N_c+1):N_c+1:(N_c+1)*(N_c+1),N_c*N_c+N_c+2:(N_c+1)*(N_c+1)-1];
A=sparse(I,J,K,nNodes,nNodes);
A(b,:)=0; A(:,b)=0; F(b)=0; A(b,b)=speye(length(b),length(b));


disp('Finish Assembly');
% solve
u = A\F;
disp('Finish Solving Au=f')
toc
Z = reshape(u, N_c+1, N_c+1)';
%[X, Y] = meshgrid(x, y);
%figure(1);
% = surf(X, Y, Z);
% title('Node value')
% set(h,'edgecolor','none');

figure(1)
% nSamples = 10;
for i = 1:N_c
    for j = 1:N_c
        % basefun(X, Y, m, n, x, y, node, type)
        vs = cellVertices(x, y, i, j);
        xlow = vs(1, 1);
        xhigh = vs(3, 1);
        ylow = vs(1, 2);
        yhigh = vs(3, 2);
        xs = linspace(xlow, xhigh, nSamples);
        ys = linspace(ylow, yhigh, nSamples);
        zs = zeros(nSamples, nSamples);
        for node = 1:4
            gloidx = loc2glo(N_c, i, j, node);
            nodevalue = u(gloidx);
            [B,a]=basefun(x, y, i, j,N_f);
            for p = 1:nSamples
                for q = 1:nSamples
                    zs(p, q) = zs(p, q) + nodevalue* a((q-1)*nSamples+p,node);
                end
            end
        end
        surf(xs, ys, zs');
        hold on;
    end
end
% title('Interpolation Field');
% memory
