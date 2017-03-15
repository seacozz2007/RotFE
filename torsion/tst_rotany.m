%tst_rotdraw.m
%
% demonstrate building a rotor model (as a Matlab structure)
% and plotting it in (fancy) 3d plot
%
% By: 	Izhak Bucher
% Date:	8-June-2000
% How: 	Free for any non-profit use (no commercial use allowed)
%      As the author of this software I specifically object to
%      commerial bodies distributing this software from their  
%      WEB and forcing users to register.
% Where: mebucher@tx.technion.ac.il


R.NODES=[0:.1:1];		% define nodal points for model
R.ELEMENTS=[];
for q=1:length(R.NODES)-1	% construct element definition
   R.ELEMENTS=[R.ELEMENTS	% each line:
      % node1 node2  d_out d_in material_no
      q q+1 (1+q*(q<5)+2*(q>8))*20e-3 16e-3 1];
end
   
%               no   E   rho   nu   <== Material definition card
   R.MATERIALS=[1   2e11 7800 .3;  2 2e11/3 7800/3 .3];
   
   %        node  d_out d_in material_no  (DISCS... or DISKS)
   R.DISCS=[3 0.2 0.01 2e-2 2	
      8 0.2 0.01 2e-2 2];
   
   R.POINT_MASS=[7 2 ];		% define a point mass at node 7, m=2(kg)
   R.SPRINGS=[];				% no springs for torsional vibration (although shaffet2 supports it, this is a give away demo so call me if you're into vibration)
   
   
   %  use the model definition to plot the model
   % note that the springs are drawn ion the direction of principal axes
   drawrot3dexp(R), title(' model for which torsional vibs. will be animated')
   
   disp(' This demo takes a long time to create')
   disp(' be patient and animated torsional modes will be shown')
   R=shaffet2(R); % build FE model for torsional vibration (R.M, R.K are created)
   tanimate(R,2);
   