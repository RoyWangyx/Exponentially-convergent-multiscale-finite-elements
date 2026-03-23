N_c=8;
N_f=8;
k0 = 2;
N=N_c*N_f;
[result,K,C] = FEM(N,k0);
L=zeros(7,1);
H=zeros(7,1);
%for i=1:7
i = 3;
erro=result-MsFEM(N_c, 3,N_f,k0);
L(i)=sqrt(erro.'*C*conj(erro))/sqrt(result.'*C*conj(result))
H(i)=sqrt(erro.'*K*conj(erro)/(result.'*K*conj(result)))
%end
%save eg3_method2_m1to7_constrast_2_to_10 L H;