function M=tanimate(Rot,no,N,T)
%M=tanimate(Rot,no,[N,T])
%
%
% animate the torional vibration of a rotor system
% using the 3d plot routine
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
   if q==0, ax=axis;
   else, axis(ax); end
      if cT, view(T); end

   M(:,q+1)=getframe;
end
movie(M,4);