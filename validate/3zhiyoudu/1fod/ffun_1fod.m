function [dz]=ffun(t,Z)
%振动方程组dJSFunZ_W_R_M_E()
%m1--k1--m2--k2--m3
%f1------f2------f3
%f1=sin(wt)  f2=0  f3=0
%c1------c2------c3

%方程:
%m1*x1..+c1*x1.+k1x1=f1

%x1..=dz(2)  x1.=dz(1)=Z(2);  x1=Z(1);


init_1fod;

dz=zeros(n*2,1);
dz(1)=Z(2);
dz(2)=(-c1*Z(2)-k1*Z(1)+sin(w*t)*f1)/m1;

end