function Kst=rotkst(L,rho,r,ri,nu)
% K1=rotstiff(L,E,r,ri,nu)
% %
%   Skew symmetric part of stiffness  (upper right corner 4 x 4)
%   Due to speed variation
%
% Total KstT=dW(t)/dt*[0 Kst ; -Kst 0];
%
%  added the timoshenko shear using Kramer's formulation
% 

%  phi = E*I/(kappa*A*G*L^2)
%      G=E/2/(1+v)
%
 
%        
 
   Kst=zeros(4);
  if r==0,  return, end
    

f=ri/r;
d=2*r;  A=pi*d^2/4*(1-f^2);
 I=pi*d^4/64*(1-f^2);  Jt=I;  Jp=2*Jt;
kappa=shearc(nu,d,2*ri);   
phi=24*(1+nu)*I/(kappa*A*L^2);

     t4 = power(1.0+phi,2.0);
      t5 = 1/t4;
      t7 = 1/L*rho*Jp*t5;
      t8 = L*L;
      t9 = t8*L;
      t10 = 1/t9;
      t11 = rho*t9;
      t12 = t11*Jp;
      t13 = Jp*phi;
      t14 = t11*t13;
      t17 = t10*(-42.0*t12+210.0*t14)*t5;
      t18 = 1/t8;
      t19 = rho*t8;
      t20 = t19*Jp;
      t21 = t19*t13;
      t24 = t18*(-42.0*t20+210.0*t21)*t5;
      t25 = phi*phi;
      t27 = t11*Jp*t25;
      t30 = t18*(-56.0*t12-140.0*t27-70.0*t14)*t5;
      t33 = t18*(-210.0*t21+42.0*t20)*t5;
      t36 = t18*(70.0*t14+14.0*t12-70.0*t27)*t5;
      t39 = t10*(-210.0*t14+42.0*t12)*t5;
      Kst(1,1) = 3.0/5.0*t7;
      Kst(1,2) = -t17/840;
      Kst(1,3) = -3.0/5.0*t7;
      Kst(1,4) = -t17/840;
      Kst(2,1) = -t24/840;
      Kst(2,2) = -t30/840;
      Kst(2,3) = -t33/840;
      Kst(2,4) = -t36/840;
      Kst(3,1) = -3.0/5.0*t7;
      Kst(3,2) = -t39/840;
      Kst(3,3) = 3.0/5.0*t7;
      Kst(3,4) = -t39/840;
      Kst(4,1) = -t24/840;
      Kst(4,2) = -t36/840;
      Kst(4,3) = -t33/840;
      Kst(4,4) = -t30/840;
