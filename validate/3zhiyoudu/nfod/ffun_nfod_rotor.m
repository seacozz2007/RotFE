function [dz]=ffun_nfod(t,Z)
%振动方程组dJSFunZ_W_R_M_E()
%m1--k1--m2--k2--m3
%f1------f2------f3
%f1=sin(wt)  f2=0  f3=0
%c1------c2------c3

%方程:
%m1*x1..+c1*x1.+k1x1=f1

%x1..=dz(2)  x1.=dz(1)=Z(2);  x1=Z(1);


init_nfod_rotor;

[na n]=size(M);
fd=zeros(n,1);
fd(1,1)=cos(w*t)*f1;
invM = inv(M);

A=Z(1:n);
B=Z((n+1):(n*2));

dz=zeros(2*n,1);
dz(1:n) = B;
dz((n+1):(n*2)) = (- invM*C*B - invM*K*A +invM*fd);

end