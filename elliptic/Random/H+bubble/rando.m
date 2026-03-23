function u=rando(N)
gm=gmdistribution(zeros(1,N^2),eye(N^2));
s=rng;
r=random(gm);
u=abs(r)+0.5;
u=reshape(u,N,N);
end