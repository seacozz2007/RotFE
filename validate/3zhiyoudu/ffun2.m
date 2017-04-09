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

n=3;

fd=[cos(w*t)*f1 0 0]';

% n=6;
% E0=zeros(n/2);
% M = [M E0;E0 M];
% K = [K E0;E0 K];
% C = [C E0;E0 C];
% fd = [cos(w*t)*f1 0 0 0 0 0]';

    invM = inv(M);
    
    A=Z(1:n);
    B=Z((n+1):(n*2));
    
    dz=zeros(2*n,1);
    dz(1:n) = B;
    dz((n+1):(n*2)) = (- invM*C*B - invM*K*A +invM*fd);

%    dz((n+1):(n*2)) = (- C*B - K*A +fd);
    
%     dz(1)=Z(4);
%     dz(2)=Z(5);
%     dz(3)=Z(6);
%     dz(4)=(-c1*Z(4)-k1*(Z(1)-Z(2))+fd(1))/m1;
%     dz(5)=(-c2*Z(5)-k1*(Z(2)-Z(1))-k2*(Z(2)-Z(3)))/m2;
%     dz(6)=(-c3*Z(6)-k2*(Z(3)-Z(2)))/m3;
end