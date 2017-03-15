function M=tanimate(Rot,no,N,T)
%M=tanimate(Rot,no,[N,T])
%
% By: 	Izhak Bucher
% Date:	8-June-2000
% How: 	Free for any non-profit use (no commercial use allowed)
%      As the author of this software I specifically object to
%      commerial bodies distributing this software from their  
%      WEB and forcing users to register.
% Where: mebucher@tx.technion.ac.il
% Purpose: animate the torional vibration of a rotor system
%  	   using the 3d plot routine
%
% Usage 

if nargin<3, N=8; end
if nargin<2, no=2; end
if nargin>3, cT=1; else, cT=0; end
Rot=shaffet2(Rot);
[V D]=eig(Rot.K,Rot.M);
for q=1:length(V)
   V(:,q)=V(:,q)/max( abs (V(:,q)) );
end

drawrot3dmode(Rot,V(:,1));
      if cT, view(T); end
M=moviein(N);
for q=0:N-1
   w=sin(q*2*pi/N);
   drawrot3dmode(Rot,V(:,no)*w);
   set(gcf,'rend','opengl')
	title(sprintf('Frame %d of %d',q+1,N))
   if q==0, ax=axis;
   else, axis(ax); end
      if cT, view(T); end

   M(:,q+1)=getframe;
end
movie(M,4);
