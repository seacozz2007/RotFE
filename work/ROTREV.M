function [Wx,Wy]=rotrev(R,r)
%
%


dp=2*pi/20;
Wx=[];Wy=[];

for p=[0:dp:2*pi-dp]
   W = real(r)*cos(p)-imag(r)*sin(p);
   [vx,vy]=rotxy(R,W);
   Wx=[Wx vx];
   Wy=[Wy vy];
   
end
 

