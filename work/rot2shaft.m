function Rs=rot2shaft(R)
%Rs=rot2shaft(R)
%
% Convert the Rot structure (the format rotfe accepts)
% into a format that is acceptable by the shaft.m program
% this allows one to re-edit the shaft elements.
%
% By:		I. Bucher
% Date:	30.8.2000
% Email:	bucher@technion.ac.il
% Rev.:	0.1
% Bla:	This program is not in the public domain
%			and cannot be used by unauthorized people 
%			without a written concent from the author

   ne=size(R.ELEMENTS,1);
    elm=[];
    for q=1:ne
       z1=R.NODES(R.ELEMENTS(q,1));
       z2=R.NODES(R.ELEMENTS(q,2));
       di=(R.ELEMENTS(q,4));
       do=(R.ELEMENTS(q,3));
       E=R.MATERIALS(R.ELEMENTS(q,5),1);
       rho=R.MATERIALS(R.ELEMENTS(q,5),2);
       nu=R.MATERIALS(R.ELEMENTS(q,5),3);
       elm=[elm;	z1 z2 do di E nu rho];
    end
    
    Rs.ELEMENTS=elm;
    
