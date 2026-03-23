k0 = 9;
N=512;
[result,K,C] = FEM(N,k0);
results=FEM(2*N,k0);
L=zeros(7,1);
H=zeros(7,1);
tmp=reshape(results,2*N+1,2*N+1);
tmp=tmp(1:2:end,1:2:end);
tmp=reshape(tmp,(N+1)^2,1);
for i=1:2
erro=result-tmp;
L(i)=sqrt(erro.'*C*conj(erro))/sqrt(result.'*C*conj(result))
H(i)=sqrt(erro.'*K*conj(erro)/(result.'*K*conj(result)))
end