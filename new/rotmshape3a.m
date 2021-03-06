function [wx,wy,z,v2]=rotmshape(v,Rot,NP,pQ,cQ,noplot)
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
%  z  - coordinate for which Wu,Wv are computed.
%
% By:		I. Bucher
% Date	15-8-98
% Modified plot format  IB 7-9-1998 
if nargin<5, cQ=1;end
if nargin<4, pQ=1; end
if nargin<3, NP=10; end
PlotFormat=2;
[Wu,Wv,z]=rotmshape(v,Rot,NP);
[wx wy]=animate2(Wu,Wv,12);
wx=0.95*wx;  wy=0.95*wy;
if nargin>5, return, end
if PlotFormat==1,
   plot3(z,wx,wy,'-y','linewidth',2);
else
   plot3(z,wx(:,1),wy(:,1),'-b','linewidth',4); hold on
   dx=max(z)-min(z);
   plot3(z,0*wx(:,1),0*wy(:,1),'-.r','linewidth',1); hold on, 
   hpQ=plot3([z z]',[1;0]*wx(:,1)',[1;0]*wy(:,1)','-m','linewidth',.7); hold on
   hpQ1=surf([z z]',[1;0]*wx(:,1)',[1;0]*wy(:,1)');
   set(hpQ1,'FaceLighting','phong','FaceColor','interp',...
      'AmbientStrength',0.5,'EdgeColor','interp')
   
             light('Position',[12 1/2 1/2],'Style','infinite');
             
end
N=length(z); p=fix(N/10); p=1; 
[mm nn]=size(wx);
hold on
hcQ=plot3( ones(nn,1)*z(1:p:N)' , wx(1:p:N,:)', wy(1:p:N,:)','color',[0 .5 0],'linewidth',0.4);
 
hold off, 
ax=axis;
 MAX=max(ax(3:6));  MIN=min(ax(3:6));  
 
   axis([ax(1:2) MIN MAX MIN MAX])
   figure(gcf)
   if ~pQ, set([hpQ(:);hpQ1(:)],'visible','off'); end
   if ~cQ, set(hcQ,'visible','off'); end