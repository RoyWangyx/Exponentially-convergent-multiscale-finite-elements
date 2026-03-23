function [e_L,e_H]= err(N_c,N_e,N_f)

N=N_c*N_f;
[result,K,C] = FEM(N);
erro=result-MsFEM(N_c, N_e,N_f);
e_L=sqrt(erro'*C*erro)/sqrt(result'*C*result);
e_H=sqrt(erro'*K*erro)/sqrt(result'*K*result);


end