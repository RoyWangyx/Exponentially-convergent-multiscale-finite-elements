function a = boundary(N_c, N_c1,N_x,N_y, m, n)
%identifying boundary; 
% N: number of elements in one direction
b=[1:N_x+1,N_x+2:N_x+1:(N_x+1)*(N_y+1),2*(N_x+1):N_x+1:(N_x+1)*(N_y+1),N_x*N_y+N_y+2:(N_x+1)*(N_y+1)-1];
a = b;
if n > 1 
    b = setdiff(b,[1:N_x+1]);
end
if n<N_c1
    b = setdiff(b,[N_x*N_y+N_y+1:(N_x+1)*(N_y+1)]);
end
if m>1
    b = setdiff(b,[1+N_x+1:N_x+1:(N_x+1)*(N_y+1)]);
end
if m<N_c
    b = setdiff(b,[(N_x+1)+N_x+1:N_x+1:(N_x+1)*(N_y+1)]);
end
a = setdiff(a,b);
end