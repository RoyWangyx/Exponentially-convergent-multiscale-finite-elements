function [R,P,bub] = restrict(X, Y, m, n,N_f,t,k0)
%compute the matrix corresponding to the oversampled patch and its
%restriction operator on the edge, based on the stiffness matrix assembled
%in basefun1, and use local subproblem to solve for the bubble and harmonic
%parts.

%P is the inner product corresponding to the oversampled norm, R is the
%restriction operator
fine=(X(2)-X(1))/N_f;
N_c=length(X)-1;

[A,N_x,N_y,x,y]=basefun1(X, Y, m, n,N_f,t,k0);
if t ==1
    b = bc(N_c-2,N_c-1, N_x,N_y, m-1, n);
else
    b = bc(N_c-1,N_c-2, N_x,N_y, m, n-1);
end
%assemble boundary b

%solve the basis function corresponding to the harmonic part, whose dof
%corresponds to the boundary values
%pay attention to the case of boundary patches
if t==1 && 2<m && m<N_c-1 && 1<n && n<N_c-1 || t==2 && 2<n && n<N_c-1 && 1<m && m<N_c-1
    f=zeros(2*(N_x+N_y),2*(N_x+N_y)-1);
    f(1:2*(N_x+N_y)-1,:)=speye(2*(N_x+N_y)-1);
else
    f=speye(length(b));
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
        xf=(x(i+1)+x(i))/2;            
        yf=(y(j+1)+y(j))/2;    
        xlow = x(i);            
        xhigh = x(i+1);
        ylow = y(j);   
        yhigh = y(j+1);
        leng = -xlow+xhigh;           
        f(1) = ffun(xf,yf)*leng^2/4+(gfun(xf,ylow,k0)+gfun(xlow,yf,k0))*leng/2;            
        f(2) = ffun(xf,yf)*leng^2/4+(gfun(xf,ylow,k0)+gfun(xhigh,yf,k0))*leng/2;            
        f(3) = ffun(xf,yf)*leng^2/4+(gfun(xf,yhigh,k0)+gfun(xhigh,yf,k0))*leng/2;           
        f(4) = ffun(xf,yf)*leng^2/4+(gfun(xf,yhigh,k0)+gfun(xlow,yf,k0))*leng/2;  
        G(loc2glo(N_x, i, j, 1))= G(loc2glo(N_x, i, j, 1))+f(1);
        G(loc2glo(N_x, i, j, 2))= G(loc2glo(N_x, i, j, 2))+f(2);
        G(loc2glo(N_x, i, j, 3))= G(loc2glo(N_x, i, j, 3))+f(3);
        G(loc2glo(N_x, i, j, 4))= G(loc2glo(N_x, i, j, 4))+f(4);
    end
end
G(b)=0; 
bub=B\G;


%compute inner product of energy
%P=fine^2*speyes(N_x*N_y+2*(N_x+N_y));
%P(2*(N_x+N_y),2*(N_x+N_y))=harm'*A*harm;
P=harm.'*A*conj(harm);
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
if norm(bub)==0
    for i=1:leng-2
        bub(i) = i*(leng-1-i)/(leng-1)^2;
    end
end
end

    


