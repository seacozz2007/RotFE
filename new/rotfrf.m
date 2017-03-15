function [frfo]=rotfrf(Rot,w)
%H=rotfrf(Rot,w)
%
% Rot - rotor structure
% w =[w1 w2 ..] [Rad/s]
% I. Bucher 
%

if nargin<2, error(' frequency vector not supplied'), end
ir = Rot.RESP_DOF;
ie=  Rot.FORCE_DOF;

METHOD=1;

if METHOD==1,
   M=full(Rot.M); G=full(Rot.G); D=full(Rot.D); K=full(Rot.K);
   W=Rot.W; N=size(M,1);
   z=zeros(size(M));
   A=[ M z ; z K ]; B=[W*G+D K ; -K z ];                      % see Lee page 156
   [r lambda]=eig(B,A);                            % 2N modes
   lm=A*r;	l=inv(lm).';                           	% see lancaster
   lambda=-diag(lambda);   % a matter of definition
   
   if reducedQ(Rot),
      nf=Rot.dim-length(Rot.Bcord);
      zl=eye(2*nf); zr=zl;
      zl=zl(ie,:)'; 
      zr=zr(ir+N,:);
      Lm=l*kron(eye(2),Rot.T.')*zl; 
      r=zr*kron(eye(2),Rot.T)*r; 
   else
      Lm=l(ie,:)';       
      r=r(ir+N,:);
    end
      frf=[];  
   for q=1:length(w),            % loop on all 'measured' frequencies
      Hw=r*diag(( 1./( i*w(q) - lambda) ))*Lm;
      frf=[frf ; Hw(:).'];
   end
else,
   error(' Must use METHOD=1')
   % Need some debuggggging here
   nM=length(Rot.M); z=zeros(nM);
   A=[ Rot.M z ; z Rot.K]; B=[Rot.W*Rot.G+Rot.D Rot.K ; -Rot.K z ];                      % see Lee page 156
   [r lambda]=eig(full(B),full(A));                   % 2N modes
   lm=A*r;	l=inv(lm).';                   % see lancaster
   %ar=diag(l.'*A*r); br=diag(l.'*B*r);
   lambda=-diag(lambda);   % a matter of definition
   
   %======= suppose we measure the mode shapes (displacement), we thus have
   %	only l1,r1   for R=[r1*diag(lambda1) conj(r1)*diag(conj(lambda1)); r1 conj(r1)]
   %	and 	      L=[l1*diag(-lambda1) conj(l1)*diag(conj(-lambda1)); l1 conj(l1)]
   
   ind=find(imag(lambda) > 0);					% take positive frequencies only
   lambda1=lambda(ind); r1=r(nM+1:2*nM,ind); l1=l(nM+1:2*nM,ind); 	% that's what we measure
   
   % create "equal" length for better conditioning
   for q=1:nM,
      r1=r1/( norm(r1(:,q)) );  
      l1=l1/( norm(l1(:,q)) );  
   end
   % reconstruct the full eignvectors right and left
   
   R=[r1*diag(lambda1) conj(r1)*diag(conj(lambda1)); r1 conj(r1)];
   L=[l1*diag(-lambda1) conj(l1)*diag(conj(-lambda1)); l1 conj(l1)];
   
   % now normalize modes according to mass and stiffness
   
   modal_mass=diag(L.'*A*R);
   for q=1:nM,
      th=-lambda(q)^2*l1(:,q).'*Rot.M*r1(:,q)+l1(:,q).'*Rot.K*r1(:,q);
      r1(:,q)=r1(:,q)/sqrt( th );
      l1(:,q)=l1(:,q)/sqrt( th );
   end
   
   %==================>  Compute FRF
   Lm=l1(ie,:);
   nr1=size(r1,1); nLm=size(Lm,1);
   frf=zeros(length(w),nr1,nLm) ;  
   for q=1:length(w),            % loop on all 'measured' frequencies
      Hw=r1*diag(( -lambda1./( i*w(q) - lambda1) ))*Lm.';
      Hw=Hw+conj(r1)*diag(( -conj(lambda1)./( i*w(q) - conj(lambda1)) ))*conj(Lm.');
      frf(q,:,:)= Hw;
   end
end

if nargout>0
   frfo=frf;
else
   semilogy(w,abs(frf(:,1)))
   figure(gcf)
end



function a=reducedQ(Rot)
%
     a=0;
      if isfield(Rot,'Reduct')
         if Rot.Reduct.flag,
            a=1;
         end
      end
      