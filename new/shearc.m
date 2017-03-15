function kappa=shearc(nu,do,di)
%kappa=shearc(nu,do,di)
%
% shear coefficient for round (hollow) shafts
%
% I. Bucher 3/97

m=di/do;

kappa=6*(1+nu)*(1+m*m)^2/( (7+6*nu)*(1+m*m)^2+(20+12*nu)*m*m);

