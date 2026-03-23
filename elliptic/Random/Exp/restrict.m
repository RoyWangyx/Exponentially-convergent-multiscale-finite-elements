function [R,P,bub] = restrict(X, Y, m, n,N_f,t)
%compute the matrix corresponding to the oversampled patch and its
%restriction operator on the edge, based on the stiffness matrix assembled
%in basefun1, and use local subproblem to solve for the bubble and harmonic
%parts.

%P is the inner product corresponding to the oversampled norm, R is the
%restriction operator
fine=(X(2)-X(1))/N_f;
N_c=length(X)-1;

[A,N_x,N_y]=basefun1(X, Y, m, n,N_f,t);
b=[1:N_x+1,N_x+2:N_x+1:(N_x+1)*(N_y+1),2*(N_x+1):N_x+1:(N_x+1)*(N_y+1),N_x*N_y+N_y+2:(N_x+1)*(N_y+1)-1];
%assemble boundary b

%solve the basis function corresponding to the harmonic part, whose dof
%corresponds to the boundary values

%pay attention to the case of boundary patches
if t==1 && n==1 && 2<m<N_c-1 || t==2 && n<3 && 1<m<N_c-1
    f=zeros(2*(N_x+N_y),N_x+2*N_y-1);
    f(N_x+2:2*(N_x+N_y),:)=speye(N_x+2*N_y-1);
elseif t==1 && n==N_c-1 && 2<m<N_c-1 || t==2 && n>N_c-2 && 1<m<N_c-1
    f=zeros(2*(N_x+N_y),N_x+2*N_y-1);
    f(1:N_x+N_y,1:N_x+N_y)=speye(N_x+N_y);
    f(N_x+N_y+2:N_x+2*N_y,N_x+N_y+1:N_x+2*N_y-1)=speye(N_y-1);
elseif t==1 && m<3 && 1<n<N_c-1 || t==2 && m==1 && 2<n<N_c-1
    f=zeros(2*(N_x+N_y),2*N_x+N_y-1);
    f(2:N_x+1,1:N_x)=speye(N_x);
    f(N_x+N_y+2:2*N_x+2*N_y,N_x+1:2*N_x+N_y-1)=speye(N_x+N_y-1);
elseif t==1 && m>N_c-2 && 1<n<N_c-1 || t==2 && m==N_c-1 && 2<n<N_c-1
    f=zeros(2*(N_x+N_y),2*N_x+N_y-1);
    f(1:N_x,1:N_x)=speye(N_x);
    f(N_x+2:N_x+N_y+1,N_x+1:N_x+N_y)=speye(N_y);
    f(N_x+2*N_y+2:2*(N_x+N_y),N_x+N_y+1:2*N_x+N_y-1)=speye(N_x-1);
elseif t==1 && n==1 && m<3 || t==2 && n<3 && m==1
    f=zeros(2*(N_x+N_y),N_x+N_y-1);
    f(N_x+N_y+2:2*(N_x+N_y),:)=speye(N_x+N_y-1);
elseif t==1 && n==1 && m>N_c-2 || t==2 && n<3 && m==N_c-1
    f=zeros(2*(N_x+N_y),N_x+N_y-1);
    f(N_x+2:N_x+N_y+1,1:N_y)=speye(N_y);
    f(N_x+2*N_y+2:2*N_x+2*N_y,N_y+1:N_x+N_y-1)=speye(N_x-1);
elseif t==1 && n==N_c-1 && m<3 || t==2 && n>N_c-2 && m==1
    f=zeros(2*(N_x+N_y),N_x+N_y-1);
    f(2:N_x+1,1:N_x)=speye(N_x);
    f(N_x+N_y+2:N_x+2*N_y,N_x+1:N_x+N_y-1)=speye(N_y-1);
elseif t==1 && n==N_c-1 && m>N_c-2 || t==2 && n>N_c-2 && m==N_c-1
    f=zeros(2*(N_x+N_y),N_x+N_y-1);
    f(1:N_x,1:N_x)=speye(N_x);
    f(N_x+2:N_x+N_y,N_x+1:N_x+N_y-1)=speye(N_y-1);
else
    f=zeros(2*(N_x+N_y),2*(N_x+N_y)-1);
    f(1:2*(N_x+N_y)-1,:)=speye(2*(N_x+N_y)-1);
end


B=A;
F=-A(:,b)*f;
B(b,:)=0; B(:,b)=0;F(b,:)=f;B(b,b)=speye(length(b));
harm=B\F;

%solve the basis function corresponding to the bubble part, whose dof
%corresponds to the right hand side values

%assemble rhs
G=zeros((N_x+1)*(N_y+1), 1);
for i = 1:N_x
    for j = 1:N_y 
        for p=1:4
            G(loc2glo(N_x, i, j, p))= G(loc2glo(N_x, i, j, p))-fine^2/4;
        end
    end
end
G(b)=0; 
bub=B\G;


%compute inner product of energy
%P=fine^2*speyes(N_x*N_y+2*(N_x+N_y));
%P(2*(N_x+N_y),2*(N_x+N_y))=harm'*A*harm;
P=harm'*A*harm;
%indentifying the edge
    leng=N_f+1;
if t==1
    if m==1
        c=[(N_x+1)*N_f+1:(N_x+1)*N_f+N_f+1];
    else 
        c=[(N_x+1)*N_f+N_f+1:(N_x+1)*N_f+2*N_f+1];
    end
else
    if n==1
        c=[N_f+1:N_x+1:(N_x+1)*N_f+N_f+1];
    else
        c=[(N_x+1)*N_f+N_f+1:N_x+1:(N_x+1)*N_f*2+N_f+1];
    end
end

%intepolate
R=harm(c,:);
R=R-linspace(1,0,leng)'*R(1,:)-linspace(0,1,leng)'*R(leng,:);
R=R(2:leng-1,:);
bub=bub(c);
bub=bub-linspace(1,0,leng)'*bub(1)-linspace(0,1,leng)'*bub(leng);
bub=bub(2:leng-1);
end

    


