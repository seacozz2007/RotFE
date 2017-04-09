% Non-linear dynamic analysis of a flexible rotor supported on porous oil journal bearings
% Nonlinear dynamic analysis of a flexible rotor supported by micropolar fluid film journal bearings
% script file which defines the model of a rotor, defines                                       
% (1) geometry                                                                                  
% (2) material                                                                                  
% (3) boundary conditions                                                                       
% (4) possible reduction of the model       
% (5) various flags which affect the run which follows     
% (6) Unbalance specification
% (7) Point force participation matrices
% (8) Point mass (linear and Angular, m & J)
%
% by:seaco 20170116

%====================================================================================           
%   (1)                                                                                         
% define nodal locations                                                                        
% NODES=[z1 z2 z3 ... ]                                                                         
% units=[Length] ([Meters])                                                                     
%                                                                                               
%====================================================================================           

 NODES=[0 0.3 0.6] ;                                                              
%====================================================================================           
%   (2)                                                                                         
% define ELEMENTS                                                                               
% elements=[node1 node2 d_out d_in material_no ;  % element #1                                  
%   ...    % element #2                                                                         
%    ... ]     % etc.                                                                           
%                                                                                               
% node1, node2 : integer indices of entries in NODES                                            
% d_out, d_in  : [length][M] outer and inner diameters of a hollow shaft                        
% material _no : integer index. defines approp. row in the MATERIAL matrix (table)              
%====================================================================================           

 ELEMENTS=[ 	1 2  12e-3 0 1                                                                  
                2 3  12e-3 0 1
  		 ];                                                                          

% or for example define uniform shaft                                                           
%  N_NODES=length(NODES);                                                                       
% for q=1:N_NODES-1, ELEMENTS(q,:)=[q q+1 15e-3 0 1]; end                                       
%====================================================================================           
%   (3)                                                                                         
% define material sets                                                                          
% MATERIALS=[E1 rho1 nu1; E2 rho2 nu2; ..... ]                                                          
% units=   E - [Pa] [psi] (young modulus)                                                       
%  rho [Kg/m^3] [lb/in^3] (density)  
%  nu [.] poison ration (0.3 for most metals)                                                           
%                                                                                               
%====================================================================================           
 
 MATERIALS=[2.067e11 7800 0.3;];     
%====================================================================================           
%   (4)                                                                                         
% define discs (rigid)                                                                          
% DISCS=[node1 d_out d_in width material_no; % define disc #1                                   
%        node2 ... ]       % etc.                                                                   
% units=   node1- integer, nodes number                                                         
%   d_out, d_in - (length) [m] [in] (diameter)                                                  
%   width   - (length) [m] [in] (disc's width)                                                  
%   material_no - integer, dfines the material set                                              
%                                                                                               
%====================================================================================           
                                                                                                
 DISCS=[];                                                                         
                                                                                                
%====================================================================================           
%   (5)                                                                                         
% define boundry conditions - springs                                                           
% SPRINGS=[node1 Kxx1 Kyy1 Kxy1 Kyx1 Ktt1 Kpp1 Ktp1 Kpt1   % bearing #1                         
%   node2 Kxx2 Kyy2 Kxy2 Kyx2 Ktt2 Kpp2 Ktp2 Kpt2   % bearing #2                                
%     ...       % etc.                                                                          
%  ]                                                                                            
%  where:                                                                                       
% linear (x,y) dofs have the following stiffness matrix                                         
%  Kxxyy=[Kxx Kxy; Kyx Kyy];                                                                    
% angular (p~=- dz/dy, t~=dz/dy)                                                                
%  Kpptt=[Kpp Kpt; Ktp Ktt]                                                                     
%                                                                                               
%  units: node - positive integer, index of node location z=NODES(node)                         
%  K- spring rate  Kg/m,  lb/in etc.                                                            
%                                                                                               
%====================================================================================           
                                                                                                
  SPRINGS=[];

%====================================================================================           
%   (6)                                                                                         
% define boundry conditions dashpots                                                            
% DASHPOTS=[node1 Cxx1 Cyy1 Cxy1 Cyx1 Ctt1 Cpp1 Ctp1 Cpt1   % bearing #1                        
%   node2 Cxx2 Cyy2 Cxy2 Cyx2 Ctt2 Cpp2 Ctp2 Cpt2   % bearing #2                                
%  %====================================================================================           
%   (6)                                                                                         
% define boundry conditions dashpots                                                            
% DASHPOTS=[node1 Cxx1 Cyy1 Cxy1 Cyx1 Ctt1 Cpp1 Ctp1 Cpt1   % bearing #1                        
%   node2 Cxx2 Cyy2 Cxy2 Cyx2 Ctt2 Cpp2 Ctp2 Cpt2   % bearing #2                                
%     ...       % etc.                                                                          
%  ]                                                                                            
%  where:                                                                                       
% linear (x,y) dofs have the following stiffness matrix                                         
%  Cxxyy=[Cxx Cxy; Cyx Cyy];                                                                    
% angular (p~=- dz/dy, t~=dz/dy)                                                                
%  Cpptt=[Cpp Cpt; Ctp Ctt]                                                                     
%                                                                                               
%  units: node - positive integer, index of node location z=NODES(node)                         
%  C- dashpot rate  kg/s,                                                                       
%                                                                                               
%====================================================================================           
                                                                                                
 DASHPOTS=[];                                                                                   

%====================================================================================           
                                                                                                
  PROP_DAMP=[];  %   1 percent

%====================================================================================           
%   (7)                                                                                         
% 
% FNodeDir=[node1.dir1 ; node2.dir2; ... nodeQ.dirQ]
%   node1=1,2, ..   dir=1,2,3,4
%    dir1=1->'xx' 3->'yy'  4->'mx' 2->'my'
%

FNodeDir=[];

%====================================================================================           
%   (8)                                                                                         
% 
% RNodeDir=[node1.dir1 ; node2.dir2; ... nodeQ.dirQ]
%   node1=1,2, ..   dir=1,2,3,4
%    dir1=1->'xx' 3->'yy'  4->'mx' 2->'my'
 

RNodeDir=[];

%====================================================================================           
%   (9)                                                                                         
% boundry conditions
%
% BCNodeDir=[node1.dir1 ; node2.dir2; ... nodeQ.dirQ]
%   node1=1,2, ..   dir=1,2,3,4
%    dir1=1->'xx' 3->'yy'  4->'mx' 2->'my'
%

%BCNodeDir=[1.1 1.2 1.3 1.4 10.1 10.2 10.3 10.4] ;
BCNodeDir=[] ;

%====================================================================================           
%   (10)                                                                                         
% unbalance specification
%
% UNBALANCE=[node1 ux1 uy1 ; node2  ux2 uy2; ... nodeQ uxQ uyQ]   node1=1,2, ..
%
%  nodei specifies the node number in the finite element model
%  uxi  - m*e (mass x displacement) unbalance in the x-direction
%  uyi  - m*e (mass y displacement) unbalance in the y-direction
 

  UNBALANCE=[];
  
%====================================================================================           
%   (11)                                                                                         
% Point mass specification
%
% POINT_MASS=[node1 m Jp Jd]
%  node1=1,2, ..  m-mass [kg] , Jp,Jd - moment of inertia [Kg-m^2]
% nd1=POINT_MASS(:,1); 	% find node
% m=POINT_MASS(:,2); 	% masses
% Jp=POINT_MASS(:,3);  Jd=POINT_MASS(:,4); 

  POINT_MASS=[];
  
  

  W = 20 * pi;
  
  
  
  
  
  