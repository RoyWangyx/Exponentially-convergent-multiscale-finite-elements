function [R,P] = restrict(X, Y, m, n,N_f,t,k0)
%compute the matrix corresponding to the oversampled patch and its
%restriction operator on the edge, based on the stiffness matrix assembled
%in basefun1, and use local subproblem to solve for the bubble and harmonic
%parts.

%P is the inner product corresponding to the oversampled norm, R is the
%restriction operator
%fine=(X(2)-X(1))/N_f;
N_c=length(X)-1;

[A,N_x,N_y]=basefun1(X, Y, m, n,N_f,t,k0);
if t ==1
    b = bc(N_c-2,N_c-1, N_x,N_y, m-1, n);
    bound = boundary(N_c-2,N_c-1, N_x,N_y, m-1, n);
else
    b = bc(N_c-1,N_c-2, N_x,N_y, m, n-1);
    bound = boundary(N_c-1,N_c-2, N_x,N_y, m, n-1);
end
[tf1,loc1] = ismember(bound,b);
%assemble boundary of dof bound and Dirichlet boundary b

%solve the basis function corresponding to the harmonic part, whose dof
%corresponds to the boundary values
%pay attention to the case of boundary patches
if t==1 && 2<m && m<N_c-1 && 1<n  || t==2 && 2<n && 1<m && m<N_c-1
    f=zeros(length(b),length(bound)-1);
    f(loc1(1:length(loc1)-1),:)=speye(length(bound)-1);
else
    f=zeros(length(b),length(bound));
    f(loc1,:)= speye(length(bound));
end


B=A;
F=-A(:,b)*f;
B(b,:)=0; B(:,b)=0;F(b,:)=f;B(b,b)=speye(length(b));
harm=B\F;

%solve the basis function corresponding to the bubble part, whose dof
%corresponds to the right hand side values
%c=setdiff(1:(N_x+1)*(N_y+1),b);
%C=A; 
%G=zeros((N_x+1)*(N_y+1),N_x*N_y);
%assemble rhs

%for i = 1:N_x
 %   for j = 1:N_y 
 %       for p=1:4
 %       G(loc2glo(N_x, i, j, p),(j-1)*N_x+i)=fine^2/4;
 %       end
 %   end
%end
%C(b,:)=0; C(:,b)=0; G(b,:)=0; C(b,b)=speye(size(b),size(b));
%bub=C\G;

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
end

    


