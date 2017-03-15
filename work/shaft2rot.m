function Rs=shaft2rot(R)
%
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
       z1=R.NODES(R.ELEMENTS(:,1));
       z2=R.NODES(R.ELEMENTS(:,2));
       E=R.MATERIALS(R.ELEMENTS(:,5));
       rho=R.MATERIALS(R.ELEMENTS(:,6));
       nu=R.MATERIALS(R.ELEMENTS(:,7));
       elm=[elm;	z1 z2 do di E nu rho];
    end
    
    Rs.ELEMENTS=elm;
    
