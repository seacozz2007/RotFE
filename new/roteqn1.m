function DOF=roteqn1(NODEDIR,Bcord,dim,NOSORT)
% roteqn1     return the equation number giver the node and the dof
%
%    DOF=roteqn1(NODEDIR,Bcord,dim[,NOSORT])
%  
%    INPUT:
%
%      NODEDIR=[node1.DIR1 ; node2.DIR2; ... nodeQ.DIRQ]  
%        NODEDIR is a vector of floating point numbers
%          where
%        fix(NODEDIR) = NODE number
%         DIR=NODEDIR-fix(NODEDIR) = dof within NODE (1..4)
%          DIR-> 1=xx 2=my 3=yy 4=mx, 0=all
%      Bcord = vector of suppressed equation numbers (after suppression due to rigid connection)
%      dim = original size of matrices, before elimination of rigidly connected dofs
%
% By  I. Bucher  4-8-1998   
% Rev. 1.1
%  modified to work with version 5.0
%  modified to work with numeric rather than string input
doff=[]; DOF=[];
if isempty(NODEDIR), 
   DOF=[]; return, 
end
NODE=fix(NODEDIR); % get nodal locations
DIR=round((NODEDIR-NODE)*10);  % 0.1 -> 1 0.2 -> 2 etc.

dof=1:dim; 
if size(Bcord)>0, dof(Bcord)=[]; end
n=length(NODEDIR) ;
if n>0,
   % first scan the '.0' == all dofs for this node
   cord=[]; for q=1:n, 
      if (DIR(q)==0), 
         cord=[cord q]; 
      end, 
   end
   NODE(cord)=[]; DIR(cord)=[];
   n=n-length(cord);
   
   for q=1:n,
      doff(q)=(NODE(q)-1)*2+ DIR(q);   % single dof to eliminate   
      if DIR(q)>2, 
         doff(q)= doff(q) + (dim/2)-2;  % if in the z-related part  
      end
      
   end %for q
   for q=1:length(doff), 
      DOF(q)=find(dof==doff(q)); 
   end
else,% n==0
   DOF=1:length(dof);
end

%  add 4 entries for each node an 'xx.0' direction was stated
if length(cord)>0,
   for p=1:length(cord),
      q=cord(p);
      cord1y=[(q-1)*2+1]:[(q-1)*2+4];  
      cord1z= cord1y + dim/2;
      dof_elim=(cord(p)-1)*2+[1:2]';               % which equations to eliminate
      DOF=[DOF(:) ;  dof_elim ];                   % eliminate y-dir   
      DOF=[DOF(:) ;  dof_elim+dim/2 ];             % eliminate z-dir   
   end
end

if nargin<5,
   DOF=sort(unique(DOF));
end
