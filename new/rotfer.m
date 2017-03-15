function [Rot]=rotfee(File)  
%  Rot=rotfer('template')
% or
%  Rot=rotfer( Rot_struct)
%
% Same as rotfe, but returnes a reduced-order model
% The input argument Reduce (Rot.Reduce) detemines how many
% modes should be included in the reduction
%
 % By  I. Bucher
% Date 21-10-1998    Rev. 3.02
% Date 27-7-99       Rev. 3.03  Modified to used the Reduct sub-field (instead of Reduce)
 
if isa(File,'struct'); % input was a Rot structure
   Rot=File;			  % better use the name Rot
   clear File
   f=fieldnames(Rot);	% get all field names
   n=size(f,1);		
   for q=1:n,			% Create the variables as in a mode file
      eval([f{q} '= Rot.' f{q} ';'])
   end
   clear Rot
else
   [File ext]=nameext(File);
   FAIL=['error([''rotfe:  Errors in script file >> '' File '' << or file no found in MATLAB path''])'];   
   eval(File,FAIL); % run script to create model 
end
save_flag=0; if nargout<1, save_flag=1; end   
   
% start FE 
Nelem=size( ELEMENTS ,1);  
dof4nod=4;% dof for nod  y, dy/dx, z, dz/dx
dof1=2; % no. of dof per one direction per node
dim=(Nelem+1)*dof4nod;   % no. of unconstrained equations (free rotor) 
M=zeros( dim , dim ); K=M; G=M; D=M; KH=M; % create zero-filled matrices
Kst=M;  % newly introduced (28-8-98) pseudo stiffnes due to angular acceleration

dof=1:dim;   
%  
W=1;  % compute for normalised speed of rotation   
%>>>> NOTE that both G and KH need to be multiplied by the real W (speed of rotation)  
%  
Cd=0;  
% Note that the hysteretic damping matix D needs to be multiplied by a constant <>1
n1=ELEMENTS(:,1); n2=ELEMENTS(:,2); L=NODES(n2)-NODES(n1); 
ro=ELEMENTS(:,3)/2; ri=ELEMENTS(:,4)/2;
% hack to overcome the ambiguity in boolean & index indexing
tMATERIALS=[MATERIALS(1,:)*0; MATERIALS];	% add dummy line so that indices start from 2
mat_no=ELEMENTS(:,5); E=tMATERIALS(mat_no+1,1); rho=tMATERIALS(mat_no+1,2); 
nu=tMATERIALS(mat_no+1,3);   

for q=1:Nelem, 
   
   % create mass matrix for the element 
   Me=rotmass(L(q),rho(q),ro(q),ri(q),nu(q));
   Ke=rotstiff(L(q),E(q),ro(q),ri(q),nu(q)); % create stiffmess matrix for the element   
   Ge=rotgyro(L(q),W,rho(q),ro(q),ri(q),nu(q));% create gyro matrix for the element
   %[De KHe]= rotdamp(L(q),W,Cd);  % this is old better disreguard for now
   
    Kste=rotkst(L(q),rho(q),ro(q),ri(q),nu(q));
   % ==> conpute the indices (equation numbers) for the y & z directions  
   cord1y=[(q-1)*dof1+1]:[(q-1)*dof1+4];  
   cord1z= cord1y + dim/2;
   M(cord1y,cord1y)=M(cord1y,cord1y)+Me;   M(cord1z,cord1z)=M(cord1z,cord1z)+Me; % add mass element to y & z
   K(cord1y,cord1y)=K(cord1y,cord1y)+Ke;   K(cord1z,cord1z)=K(cord1z,cord1z)+Ke; % add stiffness
   G(cord1y,cord1z)=G(cord1y,cord1z)-Ge;   G(cord1z,cord1y)=G(cord1z,cord1y)+Ge; % skew sym!
   Kst(cord1y,cord1z)=Kst(cord1y,cord1z)+Kste;   Kst(cord1z,cord1y)=Kst(cord1z,cord1y)-Kste; % skew sym!

   cord=[cord1y  cord1z];   
  % D(cord,cord)=D(cord,cord)+De;  %
%   KH(cord,cord)=KH(cord,cord)+KHe;  %   
end
%==========================================================
% Add rigid disc elements |
%==========================================================
if exist('DISCS'),
[nd md]=size(DISCS);  
if nd>0,  %  empty ([]) entry for discs   
   nd1=DISCS(:,1); dro=DISCS(:,2)/2; dri=DISCS(:,3)/2; % node,ro,ri  
   t=DISCS(:,4);  dmat=DISCS(:,5); drho=MATERIALS(dmat,2); % width,material,rho  
   Ai=pi*(dro.^2-dri.^2); m=Ai.*t.*drho; 
   Jp=0.5*m.*(dro.^2-dri.^2); Jd=0.5*Jp+m.*t.^2/12;  
   for q=1:nd,  % add rigidly connected disks 
      nod=nd1(q);  
      cord1y=[(nod-1)*dof1+1]:[(nod-1)*dof1+2]; % y coordinates  
      cord1z= cord1y + dim/2;% z  coordinates
      M(cord1y,cord1y)=M(cord1y,cord1y) + diag([m(q) Jd(q)]); % add to mass matrix  
      M(cord1z,cord1z)=M(cord1z,cord1z) + diag([m(q) Jd(q)]); % .. and z direction  
      G(cord1y,cord1z)=G(cord1y,cord1z) + W*diag([0 -Jp(q)]); % add to gyro-matrix   
      G(cord1z,cord1y)=G(cord1z,cord1y) + W*diag([0 Jp(q)]);% ... and z
   end % for q
end % if nd>0 
else DISCS=[];
end

if ~exist('POINT_MASS'), POINT_MASS=[]; end
[nd md]=size(POINT_MASS);  
if nd>0,  %  empty ([]) entry for discs   
   nd1=POINT_MASS(:,1); 	% find node
   m=POINT_MASS(:,2); 	% masses
   Jp=POINT_MASS(:,3);  Jd=POINT_MASS(:,4); 	% polar and diametral moments o finertia
   for q=1:nd,  % add rigidly connected disks 
      nod=nd1(q);  
      cord1y=[(nod-1)*dof1+1]:[(nod-1)*dof1+2] ; % y coordinates  
      cord1z= cord1y + dim/2; % z  coordinates
      M(cord1y,cord1y)=M(cord1y,cord1y) + diag([m(q) Jd(q)]); % add to mass matrix  
      M(cord1z,cord1z)=M(cord1z,cord1z) + diag([m(q) Jd(q)]); % .. and z direction  
      G(cord1y,cord1z)=G(cord1y,cord1z) + W*diag([0 Jp(q)]); % add to gyro-matrix   
      G(cord1z,cord1y)=G(cord1z,cord1y) + W*diag([0 -Jp(q)]);% ... and z
   end % for q
end % if nd>0 
%======================================================================
% Add springs between rotor-shaft and ground (elastic bearings)   |
%======================================================================
if ~exist('SPRINGS'), nd=0; SPRINGS=[];
else, [nd md]=size(SPRINGS);  end
if nd>0,  % no  springs ([]) empty entry, so don't
   if md<9, SPRINGS(nd,9)=0; end % extend, fill zeros
   Kb=zeros(2);  
   for q=1:nd,  % add springs 
      nod=SPRINGS(q,1); 
      % linear first 
      Kb(:)=SPRINGS(q,[2 5 4 3])  ; 
      cord1yz=(nod-1)*dof1+[1 dim/2+1]; % x & y 
      K(cord1yz,cord1yz)=K(cord1yz,cord1yz) +  Kb; % add to global stiffness
      % now angular  
      Kb(:)=SPRINGS(q,[2 5 4 3]+4); 
      cord1yz=(nod-1)*dof1+[2 dim/2+2]; % x & y 
      K(cord1yz,cord1yz)=K(cord1yz,cord1yz) +  Kb; % add to global stiffness
   end % for q
end   % if nd 
%======================================================================
% Add dashpots between rotor-shaft and ground| 
%======================================================================
if ~exist('DASHPOTS'), nd=0;   DASHPOTS=[];
else, [nd md]=size(DASHPOTS); end
if nd>0,  % no  DASHPOTS ([]) empty entry, so don't   
   if md<9, DASHPOTS(nd,9)=0; end % extend, fill zeros   
   Kb=zeros(2);  
   for q=1:nd,  % add damping 
      nod=DASHPOTS(q,1);
      % linear first 
      Kb(:)=DASHPOTS(q,[2 5 4 3]);  
      cord1yz=(nod-1)*dof1+[1 dim/2+1]; % x & y 
      D(cord1yz,cord1yz)=D(cord1yz,cord1yz) +  Kb; % add to global damping  
      % now angular  
      Kb(:)=DASHPOTS(q,[2 5 4 3]+4);       cord1yz=(nod-1)*dof1+[2 dim/2+2]; % x & y 
      D(cord1yz,cord1yz)=D(cord1yz,cord1yz) +  Kb; % add to global damping  
   end % for q
end   % if nd 
%=======================================================================   
% Apply boundry conditions(remove suppressed dofs for the equations |  
%=======================================================================   
if exist('BCNodeDir'),
   if length(BCNodeDir)>0,
      cord=roteqn1(BCNodeDir,[], dim);
      M(cord,:)=[]; M(:,cord)=[]; % remove equations for dof suppressed to zero 
      K(cord,:)=[]; K(:,cord)=[]; % from M, K  and G
      G(cord,:)=[]; G(:,cord)=[];   
      D(cord,:)=[]; D(:,cord)=[]; % if D was present .. 
      KH(cord,:)=[]; KH(:,cord)=[]; % if D was present ..   
      Kst(cord,:)=[]; Kst(:,cord)=[]; % if D was present ..   

      dof(cord)=[];
   else,
      BCNodeDir=[]; cord=[];
   end
   else, cord=[];BCNodeDir=[];
end
%if exist('PROP_DAMP'),
%    if length(PROP_DAMP)>0,
%	eval('D=addpdamp(M,K,PROP_DAMP);')
%     end
% end


%==========================================================
% unbalance specification
%==========================================================

if exist('UNBALANCE'),
   [nu mu]=size(UNBALANCE);
   if nu>0,
   Fu_cos=zeros(dim,1); Fu_sin=Fu_cos;
   nod=UNBALANCE(:,1);
   cord1x=(nod-1)*dof1+1; % x & y 
 % IB 20-8-98
  Fu_cos(cord1x)=W^2*UNBALANCE(:,2);			% x - direction
  Fu_cos(cord1x+dim/2)=W^2*UNBALANCE(:,3); % and y
  Fu_sin(cord1x)=W^2*UNBALANCE(:,3);		  % x -  sine component
  Fu_sin(cord1x+dim/2)=-W^2*UNBALANCE(:,2); % y sine
  
  Fu_sin(cord,:)=[]; 	Fu_cos(cord,:)=[]; 
  end
else,
   Fu_cos =[];  Fu_sin=[]; UNBALANCE=[];
end


Bcord=cord;  % prepare output (equations which were eliminated)

%==========================================================
% force locations
%==========================================================
if exist('FNodeDir')
   FORCE_DOF=roteqn1(FNodeDir,Bcord,dim,'nosort');   
else
   FORCE_DOF=[]; FNodeDir=[];
end

%==========================================================
% response locations
%==========================================================
if exist('RNodeDir'),
   RESP_DOF=roteqn1(RNodeDir,Bcord,dim,'nosort');
else
   RESP_DOF=[]; RNodeDir=[];
end

%==========================================================
% misalignment specification
%==========================================================
if ~exist('MISTYPE')
   MISTYPE=2;  % 1-misalignment defined after connecting the rotor & bearings
end		   	% 2- misalignment defined for bearings
if exist('MISALIGN'),
   [nu mu]=size(MISALIGN);
   if nu>0,
      % K*qs=R*f
      if length( setdiff( SPRINGS(:,1),MISALIGN(:,1) ) )>0,
         warning('Rotfe: misalignment not at bearing')
      end
         tmp=MISALIGN(:,2:3); nm=sum(tmp(:)~=0);
         nmx=find(MISALIGN(:,2)~=0);
         nmy=find(MISALIGN(:,3)~=0);
         nod=[MISALIGN(nmx,1)+.1;MISALIGN(nmy,1)+.3];
         RALGN=roteqn1(nod,Bcord,dim,'nosort');

      if    MISTYPE==1,       
         R=0*M(:,1:nm); R(RALGN,1:nm)=eye(nm);
         q=R'*(K\R); 
         Fmalgn=R*((q)\[MISALIGN(nmx,2); MISALIGN(nmy,3)]);
      else
         [q tmp]=intersect(SPRINGS(:,1),MISALIGN(:,1));
         nod1=[SPRINGS(tmp,1)+.1;SPRINGS(tmp,1)+.3];
          RALGN1=roteqn1(nod1,Bcord,dim,'nosort');
          q=0*M(:,1); t=q; q(RALGN1)=[SPRINGS(tmp,2); SPRINGS(tmp,3)];
          t(RALGN1)=[MISALIGN(nmx,2); MISALIGN(nmy,3)];
          x0=(K+diag(q))\(diag(q)*t);
         Fmalgn=diag(q)*(x0+t);
      end
   end
end

%==========================================================
% save 
%==========================================================
%>>>> Clear temporary vars
cmd='clear nmx nmy tmp RALGN R nd md BC_DIR cmd m Ai Jp Jd drho dmat drhi t ext Me Ge Ke De KHe nod ro ri L  Nelem FAIL dof4nod n1 n2 Cd E nu rho q dro dri nd1 cord1y cord1z tMATERIALS dof1 mat_no Kste cord Kb cord1yz mu MISTYPE nm nod1 f n';
eval(cmd,'disp(''err'')'); % clear & no error message
clear mat_no
%<<<<<<

M=sparse(M); K=sparse(K); G=sparse(G); D=sparse(D);
%
if nargout>0,
   cmd=['save   Temp_File.tmp '] ;
   FAIL=['error([''rotfe:  cannot save matrices in file >> '' FILE '' << ''])'];   
   eval(cmd,FAIL);   
   Rot=load('Temp_File.tmp','-mat');
   delete('Temp_File.tmp');
   Rot=rmfield(Rot,{'FAIL' 'cmd'});
end


%==========================================================
% Reduce order 
%==========================================================
if exist('Reduct') % reduce model ?
   switch Reduct.Type
   case {'modal' 1}
      Rot=reduce(Rot);
   case {'mixed' 2}
      Rot=reduces(Rot);
   case {'static' 'guyan' 3}
      Rot=reduceg(Rot);   
   end
end

   
      
      
   
    [T dd]=eiglv(K,M, Reduct.Nmodes );
    Reduct.Flag=1; % mark model as reduced
       M=T'*M*T;
       K=T'*K*T;
       G=T'*G*T;
       Kst=T'*Kst*T;
       D=T'*D*T;
 end
 


if save_flag, 
   cmd=['save ' File '.mtx '] ;
   FAIL=['error([''rotfe:  cannot save matrices in file >> '' FILE '' << ''])'];   
   eval(cmd,FAIL);   
  end   
clear save_flag
