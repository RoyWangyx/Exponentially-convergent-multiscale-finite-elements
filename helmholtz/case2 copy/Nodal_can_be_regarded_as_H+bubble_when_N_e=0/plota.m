x=[0 1];
y=[0 1];
C=zeros(1025);
for i=1:1025
    for j=1:1025
        C(i,j)=afun((j-1)/1024,(i-1)/1024);
    end
end
imagesc(x,y,C)
colorbar