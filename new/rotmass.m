function M=rotmass(L,dens,r,ri,nu)
%
%new rotor (timoshenko beam) mass matrix
% 28-8-98
%
% generated by the matlab (+symbolic) function timo_symb2.m
%
% I. bucher
%
% To be used within rotfe.m
 

%r??�?   ri ??�?

f=ri/r;
d=2*r;  A=pi*(r^2-ri^2);     %/4*(1-f^2);
I=pi*d^4/64*(1-f^4); Jt=I; 

kappa=shearc(nu,d,2*ri);  rho=dens;
phi=24*(1+nu)*I/(kappa*A*L^2) ;

    t1 = L*L;
      t2 = t1*L;
      t3 = 1/t2;
      t4 = t1*t1;
      t5 = rho*t4;
      t6 = t5*A;
      t7 = A*phi;
      t8 = t5*t7;
      t9 = phi*phi;
      t10 = A*t9;
      t11 = t5*t10;
      t12 = rho*t1;
      t13 = t12*Jt;
      t17 = power(1.0+phi,2.0);
      t18 = 1/t17;
      t19 = t3*(-312.0*t6-588.0*t8-280.0*t11-1008.0*t13)*t18;
      t20 = rho*t2;
      t21 = Jt*phi;
      t22 = t20*t21;
      t24 = rho*t4*L;
      t25 = t24*t7;
      t26 = t24*t10;
      t27 = t20*Jt;
      t28 = t24*A;
      t34 = t3*(-252.0*t8-140.0*t11-108.0*t6+1008.0*t13)*t18;
      t38 = 1/t1;
      t39 = t12*t21;
      t44 = t20*Jt*t9;
      t47 = t38*(-7.0*t26-112.0*t27-8.0*t28-280.0*t44-140.0*t22-14.0*t25)*t18;
      t53 = t38*(7.0*t26+14.0*t25+28.0*t27+6.0*t28-140.0*t44+140.0*t22)*t18;
      M(1,1) = -t19/840;
      M(1,2) = -t3*(420.0*t22-77.0*t25-35.0*t26-84.0*t27-44.0*t28)*t18/840;
      M(1,3) = -t34/840;
      M(1,4) = -t3*(420.0*t22+35.0*t26+63.0*t25-84.0*t27+26.0*t28)*t18/840;
      M(2,1) = -t38*(-44.0*t6-84.0*t13-35.0*t11-77.0*t8+420.0*t39)*t18/840;
      M(2,2) = -t47/840;
      M(2,3) = -t38*(-35.0*t11-26.0*t6-420.0*t39+84.0*t13-63.0*t8)*t18/840;
      M(2,4) = -t53/840;
      M(3,1) = -t34/840;
      M(3,2) = -t3*(-420.0*t22-63.0*t25-35.0*t26+84.0*t27-26.0*t28)*t18/840;
      M(3,3) = -t19/840;
      M(3,4) = -t3*(84.0*t27+44.0*t28-420.0*t22+35.0*t26+77.0*t25)*t18/840;
      M(4,1) = -t38*(35.0*t11-84.0*t13+26.0*t6+420.0*t39+63.0*t8)*t18/840;
      M(4,2) = -t53/840;
      M(4,3) = -t38*(-420.0*t39+77.0*t8+84.0*t13+35.0*t11+44.0*t6)*t18/840;
      M(4,4) = -t47/840;
