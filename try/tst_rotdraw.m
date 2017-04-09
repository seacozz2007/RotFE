R.NODES=[0:.1:1];
R.ELEMENTS=[];
for q=1:length(R.NODES)-1
   R.ELEMENTS=[R.ELEMENTS
      q q+1 (1+q*(q<5)+2*(q>8))*20e-3 16e-3 1];
   end
   
   R.MATERIALS=[1 2e11 7800 .3
   				2 2e11/3 7800/3 .3];

R.DISCS=[3 0.2 0.01 2e-2 2
         8 0.2 0.01 2e-2 2];

R.POINT_MASS=[7 2 ];
R.SPRINGS=[1 1e6 1e6 0 0
   11 1e6 1e6 -1e5 -1e5];
