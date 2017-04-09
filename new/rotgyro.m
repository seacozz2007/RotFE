function G=rotgyro(L,W,dens,r,ri,nu) 
%  G1=rotgyro(L,W,dens,r,ri,nu) 
%%
% By 	I. Bucher
% Date	15-3-1996
% Rev.	1.0
% For	RRA

  
  G=zeros(4);
  d=2*r;   
 Omega=W; rho=dens;

f=ri/r;
d=2*r;  A=pi*d^2/4*(1-f^2);
 I=pi*d^4/64*(1-f^4);  Jt=I;  Jp=2*Jt;
kappa=shearc(nu,d,2*ri);  rho=dens;
phi=24*(1+nu)*I/(kappa*A*L^2);
 

      t4 = power(1.0+phi,2.0);
      t5 = 1/t4;
      t7 = 1/L*rho*Jp*t5;
      t8 = L*L;
      t9 = t8*L;
      t10 = 1/t9;
      t11 = rho*t9;
      t12 = Jp*phi;
      t13 = t11*t12;
      t14 = t11*Jp;
      t17 = t10*(420.0*t13-84.0*t14)*t5;
      t18 = 1/t8;
      t19 = rho*t8;
      t20 = t19*t12;
      t21 = t19*Jp;
      t24 = t18*(420.0*t20-84.0*t21)*t5;
      t25 = phi*phi;
      t27 = t11*Jp*t25;
      t30 = t18*(-112.0*t14-140.0*t13-280.0*t27)*t5;
      t33 = t18*(-420.0*t20+84.0*t21)*t5;
      t36 = t18*(28.0*t14-140.0*t27+140.0*t13)*t5;
      t39 = t10*(84.0*t14-420.0*t13)*t5;
      G(1,1) = 6.0/5.0*t7;
      G(1,2) = -t17/840;
      G(1,3) = -6.0/5.0*t7;
      G(1,4) = -t17/840;
      G(2,1) = -t24/840;
      G(2,2) = -t30/840;
      G(2,3) = -t33/840;
      G(2,4) = -t36/840;
      G(3,1) = -6.0/5.0*t7;
      G(3,2) = -t39/840;
      G(3,3) = 6.0/5.0*t7;
      G(3,4) = -t39/840;
      G(4,1) = -t24/840;
      G(4,2) = -t36/840;
      G(4,3) = -t33/840;
      G(4,4) = -t30/840;
G=G*W;
