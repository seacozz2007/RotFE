function [Roto,M,K,F]=shaffet2(Rot,cmd1,cmd2)
%[Rot,M,K,F]=shaffet2(Rot)
%  
%
% build a global mass & stifness matrices of a round shaft beam
% undergoing torsional vibration
%
% INPUT:
%   Rot   - a vector of nodal point locations x =[x0 x1 ... ]
%
% OUTPUT:
%   Rot - A structure containing Rot.M Rot.K 
%   M,K - mass and stiffness matrices
%   F	- weighting matrix at locations such that M*q''+K*q=F*f(t)
%  
% Date:	1:November:1999
% By:	Izhak Bucher
% Ver:	1.0
  
 if nargin<1
     shaffet2 demo
     return
 end
 if isstr(Rot) % demo mode
    if nargin>1, no=cmd1; if isstr(no), no=eval(no); end,else,no=3;end
    R.NODES=[0:.1:1.1];
    R.ELEMENTS=[];
    for q=1:length(R.NODES)-1
       R.ELEMENTS=[R.ELEMENTS; q q+1 20e-3 16e-3 1];
    end
    R.MATERIALS=[ 2e11 7800 .3];
    if nargin>2,  % demo #2 add discs
       R.DISCS=[1 0.05 0.01 2e-2 1
          fix(length(R.NODES)/2) 0.05 0.01 2e-2 1
          length(R.NODES) 0.05 0.01 2e-2 1];
    end
    Rot=shaffet2(R); % <<<<<========= call FE program
    
    [V D]=eig(Rot.K,Rot.M); [d jj]=sort(diag(D)); V=V(:,jj);
    Rot.v=V; Rot.d=d;
    nx=size(V,1); x=Rot.NODES;
    plot(x(1:nx)',V(:,no),x(1:nx)',-cos((no-1)*x(1:nx)/x(nx)*pi)*max(V(:,no)),'.')
    title(sprintf(' Wn=%g  [Hz]',sqrt(d(no))/2/pi))
    if nargout>0, Roto=Rot; end
    return
 end  %<<< end demo
 
 Nelem=size(Rot.ELEMENTS,1); 
   dim=max(max(Rot.ELEMENTS(:,1:2))); % Size of M & K
   M=zeros(dim); K=M;



   L=diff(Rot.NODES);   
   
n1=Rot.ELEMENTS(:,1); n2=Rot.ELEMENTS(:,2); L=Rot.NODES(n2)-Rot.NODES(n1); 
ro=Rot.ELEMENTS(:,3)/2; ri=Rot.ELEMENTS(:,4)/2;
% hack to overcome the ambiguity in boolean & index indexing
tMATERIALS=[Rot.MATERIALS(1,:)*0; Rot.MATERIALS];	% add dummy line so that indices start from 2
mat_no=Rot.ELEMENTS(:,5); E=tMATERIALS(mat_no+1,1); rho=tMATERIALS(mat_no+1,2); 
nu=tMATERIALS(mat_no+1,3);  Ge=E/(2*(1+nu));

  
for i=1:Nelem,
	doi=2*ro(i); dii=2*ri(i); Gi=Ge(i); rhoi=rho(i);
	Jp=pi*doi^4/32-pi*dii^4/32; Li=L(i);
	Me=rho(i)*Jp*Li/6*[2 1;1 2];	% mass matrix per element
	Ke=Gi*Jp/Li*[1 -1; -1 1];	% stiffness, per element
       ind=[n1(i) n2(i)];
      M(ind,ind) = M(ind,ind) + Me;	% add to global matrix
      K(ind,ind)=K(ind,ind) + Ke;
  end   

if isfield(Rot,'DISCS');
	nd=size(Rot.DISCS,1); 
	for i=1:nd,
		%[ node_i dia_i w_i rhoi]
      nodei=Rot.DISCS(i,1);  di=Rot.DISCS(i,2); 
      wi=Rot.DISCS(i,3);  d_mat_no=Rot.DISCS(i,5);
      rhoi=tMATERIALS(d_mat_no+1,2); 
		 mi=pi*(di/2)^2*wi*rhoi; Jp=1/2*mi*(di/2)^2;
		M(nodei,nodei)=M(nodei,nodei)+Jp; % add disc's inertia
 	end
end 


if isfield(Rot,'xf')
   nf=length(Rot.xf); F=zeros(mn,nf);
	for i=1:nf,
		ind=max(find(Rot.xf(i)>=x)); % first node
		fi=[1-xf(i)/L(ind); xf(i)/L(ind)];
		F(ind:ind+1,i)=fi;
	end
 Rot.F=F;   
end	

if isfield(Rot,'flex_links')
	[n m]=size(Rot.flex_links);
	if n>0,
		kf=zeros(size(K)); % prepare
		for q=1:n,
			ind=Rot.flex_links(q,1:2); 
			k=Rot.flex_links(q,5); 
			r1=Rot.flex_links(q,3); r2=Rot.flex_links(q,4); 
			kf(ind,ind)=kf(ind,ind)+[r1^2 -r1*r2; -r1*r2 r2^2]*k;
		end
		 K=K+kf;
	end
end
 if isfield(Rot,'bc')
 if length(Rot.bc)>0,
	 M(bc,:)=[]; M(:,bc)=[]; K(bc,:)=[]; K(:,bc)=[]; F(bc,:)=[];
 end
end

Rot.M=M; Rot.K=K;
if nargout>0
   Roto=Rot;
end
