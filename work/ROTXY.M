function [vx,vy]=rotxy(R,v)
% break a response vector into xz and yz planes
%
%



    Bcord=R.Bcord; nodes=R.NODES;
   
   if length(Bcord)>0,
      N=length(v)+length(Bcord);
      v2=zeros(N,1);
      z=1:N; z(Bcord)=[];  
      v2(z)=v;
   else,
      v2=v(:);
   end
   N=length(v2);
   vx=v2(1:N/2); vy=v2(N/2+1:N);
   L=diff(nodes); L=L(:); 
   Nnodes=length(nodes); Nelem=Nnodes-1;
   
   
   
