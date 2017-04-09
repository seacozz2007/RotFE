function [C,KH]=rotdamp(L,W,Cd) 
% [C,KH]=rotdamp(L,W,Cd) 
%
% returns the internal damping matrices (C*q'+KH*q)
%
% By 	I. Bucher
% Date	15-3-1996
% Rev.	1.0
% For	RRA


   TMP=zeros(4);  C=zeros(8); KH=C;

           t1 = L^2;
      TMP(1,1) = 156;
      TMP(1,2) = 22*L;
      TMP(1,3) = 54;
      TMP(1,4) = -13*L;
      TMP(2,1) = 22*L;
      TMP(2,2) = 4*t1;
      TMP(2,3) = 13*L;
      TMP(2,4) = -3*t1;
      TMP(3,1) = 54;
      TMP(3,2) = 13*L;
      TMP(3,3) = 156;
      TMP(3,4) = -22*L;
      TMP(4,1) = -13*L;
      TMP(4,2) = -3*t1;
      TMP(4,3) = -22*L;
      TMP(4,4) = 4*t1;

KH(1:4,5:8)=TMP*(L*W*Cd/420); KH(5:8,1:4)=-TMP*(L*W*Cd/420);
C(1:4,1:4)=TMP*(L*Cd/420); C(5:8,5:8)=TMP*(L*Cd/420);
