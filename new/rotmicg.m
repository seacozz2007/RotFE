function [mo,z_cg,Ip,Id]=rotmicg(Rot,arg1) 
%[m,z_cg,Ip,Id]=rotmicg(NODES,ELEMENTS,MATERIALS,DISCS,POINT_MASS) 
%   or  
%[m,z_cg,Ip,Id]=rotmicg(TEPLATE_FILE_NAME)  
%   
% compute (for a rotor, assumed rigid): 
% (1) total mass
% (2) location of cg (centre of gravity)
% (3) polar moment of inertia (Ip)   [kg-m^2] or [lb-in^2]  
% (4) diametrical moment of inertia  [kg-m^2] or [lb-in^2]  
%   
% By  I. Bucher 
% Date 29-7-1998
% Rev. 1.1  
% Rev. 1.2 By Alan Stevens AlnStevens@aol.com
%           corrected bugs & added point-mass support
% Rev. 1.21  fixed shape of vectors for point mass


PRINT=0; if nargout<1, PRINT=1; end
if isstr(Rot); % is it a file ?   
   File=Rot; Rot=rotfe(File);
elseif isa(Rot,'struct'),
%   disp(' Rot structure provided as input')
else
   error('Rotmicg:  Input must be a file name or a Rot structure'),
end 

n1=Rot.ELEMENTS(:,1); n2=Rot.ELEMENTS(:,2);
L=Rot.NODES(n2)-Rot.NODES(n1);  % length of elements   
z=0.5*(Rot.NODES(n2)+Rot.NODES(n1));   % cg of elements assuming uniform diameter  

%  
ro=Rot.ELEMENTS(:,3)/2; ri=Rot.ELEMENTS(:,4)/2;
%  

mat_no=Rot.ELEMENTS(:,5);  
E=Rot.MATERIALS(mat_no,1); rho=Rot.MATERIALS(mat_no,2);

[nm tmp]=size(Rot.POINT_MASS);
if nm>0, % concentrated masses and discs ?
   pnd1=Rot.POINT_MASS(:,1); 	% find node
   pm=Rot.POINT_MASS(:,2); 	% masses
   pJp=Rot.POINT_MASS(:,3);  pJd=Rot.POINT_MASS(:,4); 
   pz=Rot.NODES(pnd1); % locations of point masses
else
   pm=0; pJp=0; pJd=0; pz=0;
end

 

% add DISCS' (rigid) contributions 
[nd tmp]=size(Rot.DISCS);   
if nd>0,  % DISCS<>[] (if empty, i.e. [], no discs)
   nd1=Rot.DISCS(:,1); dro=Rot.DISCS(:,2)/2; dri=Rot.DISCS(:,3)/2; % node,ro,ri  
   dh=Rot.DISCS(:,4);  dmat=Rot.DISCS(:,5); drho=Rot.MATERIALS(dmat,2); % width,material,rho  
   dz=Rot.NODES(nd1); % z-location of discs   
end

%  Id=m/4*(R^2+1/3*L^2),  Ip=m/2*R^2, Id/Ip=1/2+1/6*(L/R)^2 
% m=pi*rho*R^2*L;   
 me=pi*rho.*L(:).*(ro.^2-ri.^2); % mass of individual elements  
if nd>0,   
   dme=pi*dh.*drho.*(dro.^2-dri.^2);  
else, dme=0;   
end
m=sum(me)+sum(dme)+sum(pm); % total mass   
%   

%zcg=sum(m*z)/sum(m)
mz=me.*z(:); 
if nd>0, dmz=dme.*dz(:); else, dmz=0; end 
if nm>0, dpz=pm(:).*pz(:); else dpz=0; end

z_cg=(sum(mz)+sum(dmz)+sum(dpz))/m; 

%-------------   Section below modified by Alan Stevens August 1999    -----------
%---------------------------------------------------------------------------------

%Ipe=0.5*me.*(ro-ri).^2; if nd>0; Ipd=0.5*dme.*(dro-dri).^2; else, Ipd=0; end 
% Line above is INCORRECT. Correct line below. 
Ipe=0.5*me.*(ro.^2+ri.^2); if nd>0; Ipd=0.5*dme.*(dro.^2+dri.^2); else, Ipd=0; end

%---------------------------------------------------------------------------------
%---------------------------------------------------------------------------------

Ip=sum(Ipe)+sum(Ipd)+sum(pJp); % total polar moment of inertia, shaft + discs  + masses 

% compute diametral moment of inertia + steiner for each element
Ide=0.25*me.*(ro.^2+ri.^2+1/3*L(:).^2)+me.*(z(:).^2);  
if nd>0, Idd=0.25*dme.*(dro.^2+dri.^2+1/3*dh(:).^2)+dme.*dz(:).^2; else, Idd=0; end 

%-------------   Section below modified by Alan Stevens March 2000     -----------
%---------------------------------------------------------------------------------

% Diametral moment of inertia for point masses
pJd = pJd + pm(:).*pz(:).*pz(:);		% Moment taken about "zero".
									% Assumes input value is about point of attachment.

% add all contributions and ..  
% use Steiner to transform to cg
%Id=sum(Ide)+sum(Idd)-m*z_cg^2;

Id = sum(Ide) + sum(Idd) + sum(pJd) - m*z_cg^2; % Point mass effect added

%---------------------------------------------------------------------------------
%---------------------------------------------------------------------------------

if PRINT,
   fprintf(' \n   total mass                       :  %6.4f	[M]\n',m);
   fprintf(' \n   polar moment of inertia          :  %6.4f	[M*L^2]\n',Ip);
   fprintf(' \n   diametral moment of inertia      :  %6.4f	[M*L^2]\n',Id);
   fprintf(' \n   CENTRE  of gravity at:           :  %6.4f	[L]\n',z_cg);
   
end

if nargout>0, mo=m; end


