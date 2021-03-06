% 

w = [8.8 18.9
     33.2   11.5
	 15.5 18.77
	 6.16   16
	 28.2  11.9
	 6       18.5
	 23   12
	 6   16
	 15.44  19
	 32.5 19
	 49   11.3
	 8     17.6]*1e-3;  %M not mm
	 

 [NODES ELM]=get_elem_shaft(w);
  
   
%====================================================================================           
%   (3)                                                                                         
% define material sets                                                                          
% MATERIALS=[E1 rho1 nu1; E2 rho2 nu2; ..... ]                                                          
% units=   E - [Pa] [psi] (young modulus)                                                       
%  rho [Kg/m^3] [lb/in^3] (density)  
%  nu [.] poison ration (0.3 for most metals)                                                           
%                                                                                               
%====================================================================================           
                                                                                                
 MATERIALS=[210e9 7800 0.3;210e9 7800 0.3;  70e9 3200 0.3];     
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
                
ELEMENTS=ELM;
return
                                                                        
id1=58;

id2=30;
id3=44;
id4=90;
id5=100;
     DISCS=[%	id1 0.047 	20e-3 40e-3 2  % iron disk hub
                      % id1 (0.2-0.024) 	0.047 (16-8.4)/1000 2 % iron disk web
                       % id1 0.2 	(0.3-0.024) (16)/1000 2 % iron disk web

                        id2 0.047 	20e-3 40e-3 2  % iron disk hub
                       id2 (0.2-0.024) 	0.047 (16-8.4)/1000 2 % iron disk web
                        id2 0.2 	(0.2-0.024) (16)/1000 2 % iron disk web
                        
                        id3 0.056 	20e-3 (8.5+5+8.5)/1000 3  % iron disk hub
                       id3 0.145 0.056  5e-3 3 % iron disk web
                        
                        id4 0.056 	20e-3 (8.5+5+8.5)/1000 3  % iron disk hub
                       id4 0.145 0.056  5e-3 3 % iron disk web

                       id5 0.075 	20e-3 30e-3 2  % iron disk hub
                       id5 (0.3) 	0.075 3e-3 2 % iron disk web
                        ]; 
      
  
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
                
kc=1e6;  % coupling stiffness
 kb=10e6;
  SPRINGS=[	1 kc kc   % the missing entries, e.g. Kxy etc.                                     
                        10 kb kb
                        75 kb kb];
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

FNodeDir=[ 	 	];

%====================================================================================           
%   (8)                                                                                         
% 
% RNodeDir=[node1.dir1 ; node2.dir2; ... nodeQ.dirQ]
%   node1=1,2, ..   dir=1,2,3,4
%    dir1=1->'xx' 3->'yy'  4->'mx' 2->'my'
 

RNodeDir=[ ];

%====================================================================================           
%   (9)                                                                                         
% boundry conditions
%
% BCNodeDir=[node1.dir1 ; node2.dir2; ... nodeQ.dirQ]
%   node1=1,2, ..   dir=1,2,3,4
%    dir1=1->'xx' 3->'yy'  4->'mx' 2->'my'
%

BCNodeDir=[  ] ;

  PROP_DAMP=[];  %   1 percent

%====================================================================================           
%   (10)                                                                                         
% unbalance specification
%
% UNBALANCE=[node1 ux1 uy1 ; node2  ux2 uy2; ... nodeQ uxQ uyQ]   node1=1,2, ..
%
%  nodei specifies the node number in the finite element model
%  uxi  - m*e (mass x displacement) unbalance in the x-direction
%  uyi  - m*e (mass y displacement) unbalance in the y-direction
 

  UNBALANCE=[N_NODES 1e-3 0];
  
%====================================================================================           
%   (11)                                                                                         
% Point mass specification
%
% POINT_MASS=[node1 m Jp Jd]
%  node1=1,2, ..  m-mass [kg] , Jp,Jd - moment of inertia [Kg-m^2]
% nd1=POINT_MASS(:,1); 	% find node
% m=POINT_MASS(:,2); 	% masses
% Jp=POINT_MASS(:,3);  Jd=POINT_MASS(:,4); 

  POINT_MASS=[ ];
 