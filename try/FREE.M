% | rotor FE model template file 
% | file created >> 20-Mar-97 << 
% | Matlab script file which defines the FE model of a rotor:           
% |      (1) location of nodes                                          
% |      (2) shaft elements                                             
% |      (3) material properties                                        
% |      (4) rigid disks                                                
% |      (5) spring elements                                            
% |      (6) dashpot elements                                           
% |      (7) boundary conditions                                        
% |   [  (8) possible reduction of the model   ] not yet implemented    
% |   [  (9) various flags which affect the computation of the model ]  
% |------------------------------------------------------------------------ 
TITLE=['   ']; 
% |======================================================================== 
% |      (1)                                 
% | define nodal locations along shaft axis  
% | NODES=[z1 z2 z3 ... ]                    
% |------------------------------------------------------------------------ 
NODES = [  0:0.01:1  ]; 
 
% |======================================================================== 
% |      (2)                                                             
% | define ELEMENTS                                                      
% | elements=[node1 node2 d_out d_in material_no ;  % element #1         
% |     ...       % element #2                                           
% |     ... ]     % etc.                                                 
% |                                                                      
% | node1, node2 : integer indices of entries in NODES                   
% | d_out, d_in  : [length] outer and inner diameters of a hollow shaft  
% | material _no : defines row in the MATERIAL matrix (table)            
% |                                                                      
% |      or for example define uniform shaft                             
% | N_NODES=length(NODES);                                               
% | for q=1:N_NODES-1, ELEMENTS(q,:)=[q q+1 15e-3 0 1]; end              
% |------------------------------------------------------------------------ 
N_NODES=length(NODES);    
ELEMENTS=zeros(N_NODES-1,5);
for q=1:N_NODES-1, 
  ELEMENTS(q,:)=[q q+1 0.05 0 1 ];
 end 
 
% |======================================================================== 
% |      (3)                              
% | define table of material properties   
% | MATERIALS=[E1 rho1; E2 rho2; ..... ]  
% |  E  - Young's modulus, e.g in Pa      
% | rho - density, e.g. in Kg/m^3         
% |------------------------------------------------------------------------ 
MATERIALS = [ ...
 2.15e+011 7800  0.0 
 ]; 
 
% |======================================================================== 
% |      (4)                                                           
% | define discs (rigid)                                               
% | DISCS=[node1 d_out d_in width material_no;  % define disc #1       
% |        node2 ... ];                         % etc.                 
% |                                                                    
% |    node1       - integer, the node the disc is attached to         
% |    d_out, d_in - outer diameter, inner diameter of disc            
% |    width       - disc's width                                      
% |    material_no - integer, defines row in MATERIALS table           
% |                                                                    
% | or DISCS2=[node mass J_polar J_equatorial]  (not yet implemented)  
% |------------------------------------------------------------------------ 
DISCS = [ ...
 ]; 
 
% |======================================================================== 
% |      (5)                                                                
% |  define spring elements                                                 
% |  SPRINGS=[node1 Kxx1 Kyy1 Kxy1 Kyx1 Ktt1 Kpp1 Ktp1 Kpt1   % bearing #1  
% |           node2 Kxx2 Kyy2 Kxy2 Kyx2 Ktt2 Kpp2 Ktp2 Kpt2   % bearing #2  
% |            ...                                          ] % etc.        
% |                                                                         
% | where:                                                                  
% |     linear (x,y) dofs have the following stiffness matrix:              
% |        Kxxyy=[Kxx Kxy; Kyx Kyy];                                        
% |     angular (p='my' ~dx/dz, t='mx' ~-dy/dz) dofs have:                  
% |        Kpptt=[Kpp Kpt; Ktp Ktt]                                         
% |                                                                         
% |   node - number of node at which spring is attached                     
% |   K    - spring rate, e.g. kg/m                                         
% |                                                                         
% |  REMARKS:                                                               
% |  (a) in case any entry (Kii) is inf, a rigid connection is assumed      
% |  (b) if Kxx1 (first entry after node1) is NaN,                          
% |      all dofs for this node are rigidly connected                       
% |------------------------------------------------------------------------ 
SPRINGS = [ ... 
 ]; 
 
% |======================================================================== 
% |      (6)                                                                
% |  define dashpot damper elements                                         
% |  DASHPOTS=[node1 Cxx1 Cyy1 Cxy1 Cyx1 Ctt1 Cpp1 Ctp1 Cpt1  % bearing #1  
% |            node2 Cxx2 Cyy2 Cxy2 Cyx2 Ctt2 Cpp2 Ctp2 Cpt2  % bearing #2  
% |            ...                                          ] % etc.        
% |                                                                         
% |  where:                                                                 
% |     linear (x,y) dofs have the following damping matrix:                
% |        Cxxyy=[Cxx Cxy; Cyx Cyy];                                        
% |     angular (p='my' ~dx/dz, t='mx' ~-dy/dz) dofs have:                  
% |        Cpptt=[Cpp Cpt; Ctp Ctt]                                         
% |                                                                         
% |   node - node at which dashpot damper is attached                       
% |    C   - damper rate, e.g. in kg/s                                      
% |------------------------------------------------------------------------ 
DASHPOTS= [ ... 
 ]; 
 
% |======================================================================== 
% |      (7)                                                              
% | define boundary conditions                                            
% | BC_NODE contains a list of nodes where some boundary conditions       
% |         other than free-free occur                                    
% |  BC_NODE=[node1 node1 ... node2 ...  ];                               
% |                                                                       
% | BC_DIR contains the information of which dof is acutally clamped      
% |  BC_DIR=[dof1 dof2 ... dof1 ...];                                     
% |          possilbe dof entries are:                                    
% |           'xx'  clamped in X direction                                
% |           'yy'  clamped in Y direction                                
% |           'my'  clamped in angular direction around  Y axis ( dx/dz)  
% |           'mx'  clamped in angular direction around -X axis (-dy/dz)  
% |------------------------------------------------------------------------ 
BC_NODE = [ ...
 ]; 

BC_DIR = [ ...
  ]; 
 
% |======================================================================== 
FORCE_NODE=[]; FORCE_DIR=[]; RESP_NODE=[]; RESP_DIR=[]; 
