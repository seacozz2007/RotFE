function [Wuo,Wv,z,v2]=rotmshape(v,Rot,NP)
%[Wu,Wv,z,v]=rotmshape(v,RotorModel[,NP])
%
% This function is using the interpulation functions of the FE
% model (hermit-cubic) to regenerate the deflection curve)
% which is represented by Wu (y-direction) and Wv (z-direction) of a discrete rotor FE
%
% INPUT:
%	v	-	mode-shape (or any valid vector of displacements)
%	RotorModel	- Sturcutre containing Rotor model (see rotfe.m)
%	NP	-	(optional) no. of points in curve
%
%  OUTPUT:
%   Wu, Wv, amplitude curve in the XZ adn YZ planes
%  z  - coordinate for whcih Wu,Wv are computed.
%
% By:		I. Bucher  Date	15-8-98
%        29-7-99 treats reduced models 

% =============================	> Initialize parameters
if nargin<3,  NP=8; end
x=[0:NP-1]'/NP;		% non-dimensional variable

Bcord=Rot.Bcord; nodes=Rot.NODES;
if isfield(Rot,'Reduct')
   if Rot.Reduct.flag % reduced ?
      v=Rot.T*v; % expand
   end
end

   if length(Bcord)>0,
      N=length(v)+length(Bcord);
      v2=zeros(N,1);
      z=1:N; z(Bcord)=[];  
      v2(z)=v;
   else,
      v2=v(:);
   end
 
N=length(v2);
vy=v2(1:N/2); vz=v2(N/2+1:N);
% =============================	> Initialize Material and geom. constants
L=diff(nodes); L=L(:); 

Nnodes=length(nodes); Nelem=Nnodes-1;

% =============================	> Basic shape functions
% 	 
alfa=[polyval([2 -3 0 1],x), ...	% u
      polyval([1 -2 1 0],x), ...	% u,z
      polyval([-2 3 0 0],x), ...	% v
      polyval([1 -1 0 0],x)];		% v,z      shape functions

% =============================	>  Build global displacement shape
Wu=[]; z=[]; Wv=[];
x0=0;
n1=Rot.ELEMENTS(:,1);n2=Rot.ELEMENTS(:,2);

for q=1:Rot.dim/4-1,
   
   Li=nodes(q+1)-nodes(q);
   nod=(q-1)*4;
   phi=alfa*diag([1  Li   1 Li]); 
   u=vy((q-1)*2+1:(q-1)*2+4); 
   v=vz((q-1)*2+1:(q-1)*2+4); 
   Wu=[Wu ; phi*u ];
   Wv=[Wv ; phi*v];
   %
   z=[z ; x*Li+nodes(q)];
end

Wu=[Wu ; vy(N/2-1)];
Wv=[Wv ; vz(N/2-1)];
z=[z ; nodes( max(n2) )];

if nargout<1,
   clf
   plot(z,real([Wu Wv]))
   title(['Rotor Amplitude, Model: ' Rot.File])
else
   Wuo=Wu;
end

%
 