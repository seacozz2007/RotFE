function K=rotstiff(L,E,r,ri,nu)
% K=rotstiff(L,E,r,ri,nu)
% %
% By 	I. Bucher
% Date	30-3-1997
% Rev.	1.1
% For	FV.
%
%  added the timoshenko shear using Kramer's formulation

%  phi = E*I/(kappa*A*G*L^2)
%      G=E/2/(1+v)
%
 
%        
 
   K=zeros(4);
  if r==0,  return, end
  
f=ri/r;
d=2*r;  A=pi*d^2/4*(1-f^2);
 I=pi*d^4/64*(1-f^4);  Jt=I;  Jp=2*Jt;
kappa=shearc(nu,d,2*ri);   
phi=24*(1+nu)*I/(kappa*A*L^2);
EI=E*I;
   
 phi=24*(1+nu)*I/(kappa*A*L^2);  % 12*E*I/(kappa*G*A*L^2) and G=E/2/(1+nu)
 
      t1 = L*L;
      t3 = 1/t1/L;
      t4 = EI*phi;
      t8 = power(1.0+phi,2.0);
      t9 = 1/t8;
      t10 = t3*(-10080.0*EI-10080.0*t4)*t9;
      t11 = EI*L;
      t12 = t4*L;
      t15 = t3*(-5040.0*t11-5040.0*t12)*t9;
      t18 = t3*(10080.0*t4+10080.0*EI)*t9;
      t19 = 1/t1;
      t22 = t19*(-5040.0*t4-5040.0*EI)*t9;
      t23 = phi*phi;
      t25 = EI*t23*L;
      t28 = t19*(-3360.0*t11-840.0*t25-4200.0*t12)*t9;
      t31 = t19*(5040.0*EI+5040.0*t4)*t9;
      t34 = t19*(840.0*t25-840.0*t12-1680.0*t11)*t9;
      t37 = t3*(5040.0*t12+5040.0*t11)*t9;
      K(1,1) = -t10/840;
      K(1,2) = -t15/840;
      K(1,3) = -t18/840;
      K(1,4) = -t15/840;
      K(2,1) = -t22/840;
      K(2,2) = -t28/840;
      K(2,3) = -t31/840;
      K(2,4) = -t34/840;
      K(3,1) = -t18/840;
      K(3,2) = -t37/840;
      K(3,3) = -t10/840;
      K(3,4) = -t37/840;
      K(4,1) = -t22/840;
      K(4,2) = -t34/840;
      K(4,3) = -t31/840;
      K(4,4) = -t28/840;

