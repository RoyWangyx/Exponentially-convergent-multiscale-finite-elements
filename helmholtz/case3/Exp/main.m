N_c=32;
N_f=16;
k0 = 9;
N=N_c*N_f;
[result,K,C] = FEM(N,k0);
L=zeros(7,1);
H=zeros(7,1);
for i=1:7
erro=result-MsFEM(N_c, i,N_f,k0);
L(i)=sqrt(erro.'*C*conj(erro))/sqrt(result.'*C*conj(result))
H(i)=sqrt(erro.'*K*conj(erro)/(result.'*K*conj(result)))
end
save eg3_Ritz_os_m1to7 L H result K C;