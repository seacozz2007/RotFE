function [Fx,Fy]=FxFyFun(x,y,x_dot,y_dot)
a = AFun(x,y,x_dot,y_dot);
ff_1 = -((x-2*y_dot)^2+(y+2*x_dot)^2)^0.5/(1-x^2-y^2);
G=GFun(x,y,a);
V=VFun(x,y,a);
S=SFun(x,y,a);
Fx = ff_1 * (3*x*V-sin(a)*G-2*cos(a)*S);
Fy = ff_1 * (3*y*V+cos(a)*G-2*sin(a)*S);
end
