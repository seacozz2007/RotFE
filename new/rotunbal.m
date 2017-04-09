function [U,Wv]=rotunbal(Rot,W,N)
%[U,Wv]=rotunbal(Rot,W,N)
%
%  Compute the response to unbalance at a range of speeds
%
% The co-ordinates for which Rot.RESP_DOF 
%        (see Rot.RNodeDir in themodel file) is selected
% The unbalance is defined in Rot.UNBALANCE
%
% Input:
%  Rot:  a model file or a Rot structure 
%  W  : optional, speed range in RPM W=Wmax, W=[W1 W2]
%  N:   number of speeds between 0-Wm or W1-W2
%
% By:
%  I.Bucher 24-10-98
 
 if isstr(Rot), Rot=rotfe(Rot); end  % template to Rot
 
  nM=size(Rot.M,1); z=sparse(nM,nM);
  
  nr=length(  Rot.RESP_DOF );
  
  % construct Speed vector
  if nargin<3, N=100; end
  if nargin<2, W=10000; end
     
    
  if length(W)>1, W1=W(1)*2*pi/60; W2=W(2)*2*pi/60;
   else, W1=0; W2=W(1)*2*pi/60; end
   Wv=linspace(W1,W2,N);
    U=zeros(N, nr);
    for q=1:N,
       w=Wv(q);
       a= ( -w^2*Rot.M+i*w*(w*Rot.G+Rot.D)+Rot.K);
       b=w^2*(Rot.Fu_cos-i*Rot.Fu_sin);
       if isfield(Rot,'Reduct')
          if Rot.Reduct.flag, c=Rot.T*(a\b); else, c=a\b;end
          else, c=a\b;
       end
       U(q,:)=c( Rot.RESP_DOF ).';
     end
    
   
 