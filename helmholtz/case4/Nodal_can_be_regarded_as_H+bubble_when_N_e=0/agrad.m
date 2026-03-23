function a=agrad(t,s)
%agrad(t,s)=\frac{\partial a(t,s)}{\partial x}=\frac{\partial a(s,t)}{\partial y}
P = 1.5;
epsilon = 0.05;
    a=-2*pi/epsilon*P*cos(2*pi*t/epsilon)/(2 + P*sin(2*pi*t/epsilon))^2/(2 + P*sin(2*pi*s/epsilon));
end