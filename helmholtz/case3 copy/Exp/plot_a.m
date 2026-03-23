x=0:1e-3:1; lx=length(x);
y=0:1e-3:1; ly=length(y);
z=zeros(lx,ly);

for iter_x=1:length(x)
    for iter_y=1:length(y)
        z(iter_x,iter_y)=afun(x(iter_x),y(iter_y));
    end
end

axesfontsize=16;
axeslinewidth=1.8;
linelinewidth=1.8;
patchlinewidth=1.5;
set(0,'defaultaxesfontsize',axesfontsize,'defaultaxeslinewidth',axeslinewidth,...
    'defaultlinelinewidth',linelinewidth,'defaultpatchlinewidth',patchlinewidth)

contourf(x,y,log(z)); 
colorbar('TickLabelInterpreter', 'latex');
xlabel('$x_1$','interpreter','latex');
ylabel('$x_2$','interpreter','latex');
set(gca,'TickLabelInterpreter', 'latex');
title('$\log_{10} a(x)$','interpreter', 'latex');
h=gcf;
myprint('figure_log_a(x)',h)
