function [M,K,F]=rottors(File)
%
% process a rotfe input template file
% and create a torsional vibration model
%
% I. Bucher 20.7.98

if isstr(File); % is it a file ?  
   [File ext]=nameext(File);
   save_flag=0; if nargout<1, save_flag=1; end   
   FAIL=['error([''rottors:  script file >> '' File '' << not found in MATLAB path''])'];   
   eval(File,FAIL); % run script to create model 
elseif  isa(File,'struct'); % input was a Rot structure
   Rot=File;			  % better use the name Rot
   clear File
   f=fieldnames(Rot);	% get all field names
   n=size(f,1);		
   for q=1:n,			% Create the variables as in a model file
      eval([f{q} '= Rot.' f{q} ';'])
   end
   clear Rot

else
   error(['rottors:  Input argument File must be a string ']);
end

%>>>>>> diameters
d=ELEMENTS(:,3:4);  %  diameters external & internal

%>>>>>>  material
tMATERIALS=[MATERIALS(1,:)*0; MATERIALS];	% add dummy line so that indices start from 2
mat_no=ELEMENTS(:,5); E=tMATERIALS(mat_no+1,1); rho=tMATERIALS(mat_no+1,2); 
nu=tMATERIALS(mat_no+1,3);   

MAT=[E./(2*(1+nu)) rho]; % material input for shaffet

%>>>>>>  DISCS
[nd md]=size(DISCS);  
if nd>0,  %  empty ([]) entry for discs  
   nd1=DISCS(:,1); do=DISCS(:,2); di=DISCS(:,3); % node,ro,ri  
   t=DISCS(:,4);  dmat=DISCS(:,5); drho=MATERIALS(dmat,2); % width,material,rho  
   
   DISCStor=[nd1 do t drho];  % must add di later <<<<<<!!!! 
   
else
   DISCStor=[];
end
if ~exist('TorsionBC'),  TorsionBC=[];  end


if exist('POINT_MASS')
    [nd md]=size(POINT_MASS);  
    if nd>0,  %  empty ([]) entry for discs   
        Jp=POINT_MASS(:,[1 3]); 	% find node
    end % if nd>0 
else
    Jp=[];
    
end

[M,K,F]=shaffet(NODES,d,MAT,TorsionBC,DISCStor,eps+NODES(1),[],Jp);

if nargout==1
   save temp.mtx
   M=load('temp.mtx','-mat');
   delete temp.mtx
end
