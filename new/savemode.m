function savemode(v,Rot)
%savemode(v,Rot)
%
% print the modeshape v of a rotor described by the structure Rot.
% this function

Bcord=Rot.Bcord; nodes=Rot.NODES;

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
for q=1:Nelem,
   fprintf(' Node-%3.0f  @z=%4.3f, x=%4.3g   y=%4.3g   dx/dz=%4.3g   dy/dz=%4.3g \n', ...
      q,nodes(q),vx((q-1)*2+1),vy((q-1)*2+1),vx(q*2),vy(q*2) )
   
end


