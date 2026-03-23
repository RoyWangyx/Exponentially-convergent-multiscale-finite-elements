function [result] = MsFEM(N_c, N_e,N_f,k0)
tic
% setting parameters
%N_c= 16;
%coarse mesh
%N_e=2;
%number of elements per edge
%N_f=8;
%number of fine mesh per coarse patch
nSamples = N_f+1;
fprintf('Coarse Grid %d X %d \n', N_c, N_c);
fprintf('Fine Grid %d X %d \n', N_c*N_f, N_c*N_f);

x = linspace(0, 1, N_c+1);
y = x;
N_p=N_e+1;
nNodes = (N_c+1)*(N_c+1)+2*N_c*(N_c+1)*N_p+1;

%sparse assembling
numb=(5+4*N_p)^2*(N_c)^2;
I=zeros(numb,1);J=I;K=I;
F=zeros(nNodes, 1);
count=5+4*N_p;
val=zeros(nSamples^2,count,N_c^2);

% there might be inner patches, patches on the edge, and nodal patches
for i = 1:N_c
    for j = 1:N_c  
        [value,k, f] = elementstiff(x, y, i, j,N_f, N_e,k0);
        val(:,:,(i-1)*N_c+j)=value;
        for p = 1:count
            if p<5
                global_p = loc2glo(N_c, i, j, p);
            elseif p<5+N_p
                global_p = (i-1+(j-1)*N_c)*N_p+(N_c+1)^2+p-4;
            elseif p<5+2*N_p
                global_p = (i-1+(j)*N_c)*N_p+(N_c+1)^2+p-4-N_p;
            elseif p<5+3*N_p
                global_p = (i-1+(j-1)*(N_c+1)+N_c*(N_c+1))*N_p+(N_c+1)^2+p-4-2*N_p;
            elseif p<5+4*N_p
                global_p = (i+(j-1)*(N_c+1)+N_c*(N_c+1))*N_p+(N_c+1)^2+p-4-3*N_p;
            else
                global_p = (N_c+1)*(N_c+1)+2*N_c*(N_c+1)*N_p+1;
            end
            for q = 1:count
                if q<5
                    global_q = loc2glo(N_c, i, j, q);
                elseif q<5+N_p
                    global_q = (i-1+(j-1)*N_c)*N_p+(N_c+1)^2+q-4;
                elseif q<5+2*N_p
                    global_q = (i-1+(j)*N_c)*N_p+(N_c+1)^2+q-4-N_p;
                elseif q<5+3*N_p
                    global_q = (i-1+(j-1)*(N_c+1)+N_c*(N_c+1))*N_p+(N_c+1)^2+q-4-2*N_p;
                elseif q<5+4*N_p
                    global_q = (i+(j-1)*(N_c+1)+N_c*(N_c+1))*N_p+(N_c+1)^2+q-4-3*N_p;
                else
                    global_q = (N_c+1)*(N_c+1)+2*N_c*(N_c+1)*N_p+1;
                end
                index=(5+4*N_p)^2*(N_c)*(i-1)+(5+4*N_p)^2*(j-1)+(5+4*N_p)*(p-1)+q;
                I(index)=global_p;
                J(index)=global_q;
                K(index)=k(p, q);
            end
            F(global_p) = F(global_p) + f(p);
        end
    end
end



b=zeros(1,3+N_c+4*N_c*N_p);
b(1:3+N_c+2*N_c*N_p)=[1:N_c+1,(N_c+1)*N_c+1,(N_c+1)*(N_c+1),...
    (N_c+1)^2+1:(N_c+1)^2+N_p*N_c,...
    (N_c+1)^2+1+N_c^2*N_p:(N_c+1)^2+N_p*N_c+N_c^2*N_p];
for i=1:N_p
    b(3+N_c+2*N_c*N_p+1+(i-1)*2*N_c:3+N_c+2*N_c*N_p+(i)*2*N_c)=...
        [i+N_c*(N_c+1)*N_p+(N_c+1)^2:N_p*(N_c+1):i+N_p*(N_c+1)*(N_c-1)+N_c*(N_c+1)*N_p+(N_c+1)^2,...
        i+N_c*N_p+N_c*(N_c+1)*N_p+(N_c+1)^2:N_p*(N_c+1):i+N_c*N_p+N_p*(N_c+1)*(N_c-1)+N_c*(N_c+1)*N_p+(N_c+1)^2];
end
A=sparse(I,J,K,nNodes,nNodes);
A(b,:)=0; A(:,b)=0; F(b)=0; A(b,b)=speye(length(b),length(b));


disp('Finish Assembly');
u = A.'\F;
disp('Finish Solving Au=f')
toc
result=zeros((N_c*N_f+1));

figure(1)
for i = 1:N_c
    for j = 1:N_c
        vs = cellVertices(x, y, i, j);
        xlow = vs(1, 1);
        xhigh = vs(3, 1);
        ylow = vs(1, 2);
        yhigh = vs(3, 2);
        xs = linspace(xlow, xhigh, nSamples);
        ys = linspace(ylow, yhigh, nSamples);
        zs = zeros(nSamples, nSamples);
        value=val(:,:,(i-1)*N_c+j);
        for p = 1:count
            if p<5
                global_p = loc2glo(N_c, i, j, p);
            elseif p<5+N_p
                global_p = (i-1+(j-1)*N_c)*N_p+(N_c+1)^2+p-4;
            elseif p<5+2*N_p
                global_p = (i-1+(j)*N_c)*N_p+(N_c+1)^2+p-4-N_p;
            elseif p<5+3*N_p
                global_p = (i-1+(j-1)*(N_c+1)+N_c*(N_c+1))*N_p+(N_c+1)^2+p-4-2*N_p;
            elseif p<5+4*N_p
                global_p = (i+(j-1)*(N_c+1)+N_c*(N_c+1))*N_p+(N_c+1)^2+p-4-3*N_p;
            else
                global_p = (N_c+1)*(N_c+1)+2*N_c*(N_c+1)*N_p+1;
            end
            nodevalue = u(global_p);
            valuep=reshape(value(:,p),nSamples,nSamples);
            zs = zs + nodevalue* valuep;
        end
        %bubble part
        %zs=zs+bubble(x,y,i,j,N_f,k0);
        result((i-1)*N_f+1:(i-1)*N_f+nSamples,(j-1)*N_f+1:(j-1)*N_f+nSamples)=zs;
        %surf(xs, ys, zs');
        %hold on;
    end
end
result=reshape(result,(N_c*N_f+1)^2,1);
end
