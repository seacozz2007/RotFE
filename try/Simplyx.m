
global nn
TITLE=[' simplex  ']; 
 
%>>>>>>> nn must be specified elsewhere
NODES = [ 0:nn ]/nn;
 N_NODES=length(NODES);
for q=1:N_NODES-1, 
  ELEMENTS(q,:)=[q q+1 5e-2 0 1]; 
end

MATERIALS = [ ...
 2.15e+011 7800  0.3 
 ]; 
 
 DISCS = [ ...
 ]; 
 
 SPRINGS = [ ... 
 ]; 
 
 DASHPOTS= [ ... 
 ]; 
 
 
BCNodeDir=[1.1 1.3 nn+1+.1 nn+1+.3];
 
  
% |======================================================================== 
FORCE_NODE=[]; FORCE_DIR=[]; RESP_NODE=[]; RESP_DIR=[]; 
