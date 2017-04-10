function [dz]=ffun(t,Z)
%振动方程组dJSFunZ_W_R_M_E()
%m1--k1--m2--k2--m3
%f1------f2------f3
%f1=sin(wt)  f2=0  f3=0
%c1------c2------c3

%方程:
%m1*x1..+c1*x1.+k1(x1-x2)=f1
%m2*x2..+c2*x2.+k1(x2-x1)+k2(x2-x3)=f2
%m3*x3..+c3*x3.+k2(x3-x2)=f3

%x1..=dz(2)  x1.=dz(1)=Z(2);  x1=Z(1);
%x2..=dz(4)  x2.=dz(3)=Z(4);  x2=Z(3);
%x3..=dz(6)  x2.=dz(5)=Z(6);  x2=Z(5);

init;

dz=zeros(6,1);
dz(1)=Z(2);
dz(2)=(-c1*Z(2)-k1*(Z(1)-Z(3))+cos(w*t)*f1)/m1;
%dz(2)=(-c1*Z(2)-k1*Z(1)+sin(w*t)*f1)/m1;


dz(3)=Z(4);
dz(4)=(-c2*Z(4)-k1*(Z(3)-Z(1))-k2*(Z(3)-Z(5)))/m2;

dz(5)=Z(6);
dz(6)=(-c3*Z(6)-k2*(Z(5)-Z(3)))/m3;
end