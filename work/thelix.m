function [psi,ang_at_node,rx,ry]=thelix(R,v)
%  [psi]=thelix(R,r)
%
% Transform a response vector r into helical co-ordinatr system
% such that one of the axes coincides with the minimal stiffness
% direction of each bearing.
%
% Assume that the rotor is undaped and that the co-ordinates system
% is as described in C. W. Lee
%


if isfield(R,'SPRINGS'),
   SPRINGS=R.SPRINGS;  [nd md]=size(SPRINGS); 
else
   SPRINGS=[]; nd=0;
end

if nd>0,  % no  springs ([]) empty entry, so don't
   if md<9, SPRINGS(nd,9)=0; end % extend, fill zeros
   Kb=zeros(2);  min_node=1; ang_node=[];
   for q=1:nd,  % add springs 
      nod=SPRINGS(q,1); 
      % linear first 
      Kb(:)=SPRINGS(q,[2 5 4 3])  ; 
      [vv dd1]=eig(Kb); % find principle axes
      [dd jj]=sort(abs(diag(dd1))); vv=vv(:,jj); % min stiffness
      psi.angle{q}=atan2(-vv(1,1),-vv(2,1));
      psi.node{q}=nod; if nod>min_node, min_node=nod; end, ang_node=[ang_node;nod];
      psi.stiff{q}=sort( diag(dd1)  )';      
   end % for q
end   % if nd 

if nargin>1, % r provided
   
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
   vy=v2(1:N/2); vz=v2(N/2+1:N);
   L=diff(nodes); L=L(:); 
   Nnodes=length(nodes); Nelem=Nnodes-1;
   
   na=length(psi.angle);
   
   
   b0=psi.angle{1};  % This is the initial angle
%   b1=psi.angle{2};
   ang_at_node=ones(size(nodes))*b0;
   for q=1:na-1, 
      i1=ang_node(q);  i2=ang_node(q+1); 
      b0=psi.angle{q};  % This is the initial angle
      b1=psi.angle{q+1};

      ang_at_node(i1:i2)=([i1:i2]-i1)*(b1-b0)/(i2-i1)+b0;
    end
    
    if i2<Nnodes
         ang_at_node(i2:Nnodes)=ones(size(ang_at_node(i2:Nnodes)))* ang_at_node(i2);
      end
      
      % ====> Transform r into helical co-ordinates
      rx=zeros(Nnodes,1); ry=zeros(Nnodes,1);
      vz=-i*vz;
      for q=1:Nnodes
         a=-ang_at_node(q); T=[cos(a) sin(a); -sin(a) cos(a)];
         x=vy(2*q-1);    y=vz(2*q-1);
         xya=T*[x;y];  rx(q)=xya(1); ry(q)=xya(2);
      end
      
end
