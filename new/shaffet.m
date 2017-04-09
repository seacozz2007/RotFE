function [M,K,F]=shaftfet(x,d,MAT,bc,DISCS,xf,flex_links,JPx)
%[M,K,F]=shaftfet(x,d,MAT,bc[,DISCS],[xf],[flex_links])
%  
%
% build a global mass & stifness matrices of a round shaft beam
% undergoing torsional vibration
%
% INPUT:
%   x   - a vector of nodal point locations x =[x0 x1 ... ]
%   d   - a vector of diameters (one for each element)
%   MATERIAL =[G1 rho1; G2 rho2; ..]  Gi-shear modulus [Pa]  rho - density [M/L^3]
%   bc - boundry conditions (nodes whos numbers appear in bc are suppressed to zero)
%   DISCS = [ node_i dia_i w_i rhoi; ... ; node diameter width density]
%   xf  - [x1; x2 ..]  x-locations at which a point moment is applied
%   flex_links = [node_i node_j r1 r2 k; node_i2  ...] flexible links 
%	(stiffness k) between node_i and node_j node_i has radius r1 and 
%	node_j has radius r2
%
% OUTPUT:
%   M,K - mass and stiffness matrices
%   F	- weighting matrix at locations such that M*q''+K*q=F*f(t)
%  
% Date:	14:March:1998
% By:	Izhak Bucher
% Ver:	0.1
%
% EXAMPLE #1:
%	
% x=[0:.02:1]';  d=20e-3; MAT=[8.0769e+010 7800]; 
% [M,K]=shaffet(x,d,MAT,[]);
% [V D]=eig(K,M); [d jj]=sort(diag(D)); V=V(:,jj);
% no=3; nx=size(V,1);
% plot(x(1:nx)',V(:,no),x(1:nx)',cos((no-1)*x(1:nx)/x(nx)*pi)*max(V(:,no)),'.')
% title(sprintf(' Wn=%g  [Hz]',sqrt(d(no))/2/pi))
% 
% Example #2
%%[ node_i dia_i w_i rhoi]
% x=[0:.02:1]';  d=20e-3; MAT=[8.0769e+010 7800];
% DISCS=[23 0.1 0.02 7800; 42 0.1 0.02 7800 ];
% bc=[];
% [M,K]=shaffet(x,d,MAT,bc,DISCS);
% [V D]=eig(K,M); [d jj]=sort(diag(D)); V=real(V(:,jj));
% no=2; nx=size(V,1);
% plot(x(1:nx)',V(:,no),x(1:nx)',cos((no-1)*x(1:nx)/x(nx)*pi)*max(V(:,no)),'.')
%
% EXAMPLE #3:
%	
% x=[0:.05:1]';  d=20e-3; MAT=[8.0769e+010 7800]; 
% bc=[]; xf=[0 .86]'; DISCS=[];
% [M,K,F]=shaffet(x,d,MAT,bc,DISCS,xf);
% wa=2*pi*670;  wb=2*pi*3200;
% X=(K-wa^2*M)\F(:,1)+ (K-wb^2*M)\F(:,2); 
% nx=length(X);
% plot(x(1:nx)',X)
% title(sprintf(' wa=%g  [Hz]',wa/2/pi))
 

   mn=length(x);    % no. of nodes
	n=mn-1;		% no. of elements
	M=zeros(mn); K=M;

   if size(MAT,1)==1, 
	MAT=MAT(ones(n,1),:); 	% duplicate
  end
  if size(d,1)==1, 
	d=d(ones(n,1),:); 	% duplicate
  end

   L=diff(x);    p=1; 

 for i=1:n,
	di=d(i); Gi=MAT(i,1); rhoi=MAT(i,2);
	Jp=pi*di^4/32; Li=L(i);

	Me=rhoi*Jp*Li/6*[2 1;1 2];	% mass matrix per element
	Ke=Gi*Jp/Li*[1 -1; -1 1];	% stiffness, per element

      M(p:p+1,p:p+1) = M(p:p+1,p:p+1) + Me;	% add to global matrix
      K(p:p+1,p:p+1)=K(p:p+1,p:p+1) + Ke;
      p=p+1;
  end   

 if nargin>4,
	nd=size(DISCS,1); 
	for i=1:nd,
		%[ node_i dia_i w_i rhoi]
		nodei=DISCS(i,1);  di=DISCS(i,2); wi=DISCS(i,3);rhoi=DISCS(i,4);
		 mi=pi*(di/2)^2*wi*rhoi; Jp=1/2*mi*(di/2)^2;
		M(nodei,nodei)=M(nodei,nodei)+Jp; % add disc's inertia
 	end
 end	
 
 
 % concentrated mass (polar mass moment of inertia) 
 if nargin>7
     [nd md]=size(JPx);  
     if nd>0,  %  empty ([]) entry for discs   
         nd1=JPx(:,1); 	% find node
         Jp=JPx(:,2);
         Jpd=M(:,1)*0; Jpd(nd1)=Jp;
         M=M+diag(Jpd);
     end % if nd>0 
 end
    
 

 if nargin>5,
	nf=length(xf); F=zeros(mn,nf);
	for i=1:nf,
		ind=max(find(xf(i)>=x)); % first node
		fi=[1-xf(i)/L(ind); xf(i)/L(ind)];
		F(ind:ind+1,i)=fi;
	end
 end	

 if nargin>6, % flexible links present 
	[n m]=size(flex_links);
	if n>0,
		kf=zeros(size(K)); % prepare
		for q=1:n,
			ind=flex_links(q,1:2); 
			k=flex_links(q,5); 
			r1=flex_links(q,3); r2=flex_links(q,4); 
			kf(ind,ind)=kf(ind,ind)+[r1^2 -r1*r2; -r1*r2 r2^2]*k;
		end
		 K=K+kf;
	end
end
 if length(bc)>0,
	 M(bc,:)=[]; M(:,bc)=[]; K(bc,:)=[]; K(:,bc)=[]; F(bc,:)=[];
 end
 
