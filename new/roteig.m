function [r1,lambda1]=roteig(Rot,crit)
%[r1,lambda1]=roteig( RotorModel );
%
% Compute the eigevectors and eigenvalues
% of a Rotor system,
% Input:
%   RotorModel : a structure (or an object of type rotor)
%                 which was created by one of the two commands:
%						(1)  RotorModel=rotfe('model');
%						(2)  RotorModel=load('model.mtx','-mat');
%    Whithin the model the following fields are being used
%      M,K,D,G,W (i.e. RotorModel.M, RotorModel.K ...)
%
% By:		I. Bucher
% Date	25-03-2000,
%  added an optional field (Rot.scale_eigenvectors) allowing eigenvectro scaling
% Example:
%			R=rotfe('pmet1'); R.W=1000; [eigenvectors,eigenvalues]=roteig(	R );
nM=length(Rot.M); 
if nargin>1, % compute critical speeds & modes
   [r1 lambda1]=eig(full(Rot.K),full(Rot.M)-j*full(Rot.G));
   lambda1=sqrt(diag(lambda1));
   [ lambda1 ind]=sort(abs(lambda1));
   lambda1=j*lambda1;
   r1=r1(:,ind);
   
else
   if  Rot.W<1e-4,   % canceled the usage of eigl Rot.W==0, temporary
      eval('[r1 lambda1]=eigl(Rot.K,Rot.M,min(nM,30),[1 .01 .01],30);')
      lambda1=i*sqrt(lambda1); vecRigid=[]; nvec=0;
   else
      z=zeros(nM);
       
      A=[ Rot.M z ; z Rot.K]; B=[Rot.W*Rot.G+Rot.D Rot.K ; -Rot.K z ]; % see Lee page 156
      [r lambda]=eig(full(B),full(A));                   % 2N modes
      lambda=-diag(lambda);   % a matter of definition
      %It seems a bit difficult to remove the redundent rigid body modes
      % which are created by the double-size (1st order) eignvalue problem being solved
      %therfore an accurate SVD based range space method is being used, IB 7-99
       indRigid=find(abs(lambda)<1e-2); % find rigid body modes
       if length(indRigid)>0
          [vecRigid S V ]=svd(r(nM+1:2*nM,indRigid)); % find range space of rigid modes
          nvec=sum(diag(S)/S(1,1)>1e-3);
          vecRigid=vecRigid(:,1:nvec);
       else
          nvec=0; vecRigid=[];
       end
       %
      ind=find( (imag(lambda) > 0) & (abs(lambda)>1e-2) );	% take positive frequencies only        
      																      % but not rigid body ones
      lambda1=lambda(ind); 
      r1=r(nM+1:2*nM,ind);  	% that's what we measure
   end
   [dd jj]=sort(abs(lambda1));
   r1=[vecRigid r1(:,jj)];
   if isfield(Rot,'scale_vectors')
      switch lower(Rot.scale_vectors)
      case 'unit',  for qq=1:size(r1,2), r1(:,qq)=r1(:,qq)/norm(r1(:,qq)); end
      case 'mass',  
         for qq=1:size(r1,2), 
            mr=r1(:,qq).'*Rot.M*r1(:,qq); 
            r1(:,qq)=r1(:,qq)/sqrt(mr); 
         end
      end
   else  % perform unit scaling anyway
       for qq=1:size(r1,2), r1(:,qq)=r1(:,qq)/norm(r1(:,qq)); end
   end
   
   lambda1=[0*ones(nvec,1);lambda1(jj)];
end
 