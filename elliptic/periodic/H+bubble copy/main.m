N_c=32;
N_f=32;
N=N_c*N_f;
[result,K,C] = FEM(N);
L=zeros(7,1);
H=zeros(7,1);
for i=1:7
    
erro=result-MsFEM(N_c, i,N_f);
L(i)=sqrt(erro'*C*erro)/sqrt(result'*C*result);
H(i)=sqrt(erro'*K*erro/sqrt(result'*K*result));
disp(i); disp(L(i)); disp(H(i));
end
save eg1_method2_m1to7_addbb L H;