N_c=8;
N_f=8;
k0 = 2;
N=N_c*N_f;
global goo;
global gov;
global gob;
goo=importdata('u.mat');
gov=importdata('v.mat');
gob=importdata('beta.mat');
[result,K,C] = FEM(N,k0);
L=zeros(7,1);
H=zeros(7,1);
%for i=1:7
i = 2;
erro=result-MsFEM(N_c, i,N_f,k0);
L(i)=sqrt(erro.'*C*conj(erro))/sqrt(result.'*C*conj(result))
H(i)=sqrt(erro.'*K*conj(erro)/(result.'*K*conj(result)))
%end
%save eg3_method2_m1to7_constrast_2_to_10 L H;