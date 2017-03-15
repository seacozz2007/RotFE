L=0.6;D=2e-3;E=2.067e11;
I=pi*D^4/64;
F=1;
ansY=F*L^3/3/E/I
K=F/ansY
%build the fe model
Rot=rotfe('ROTDEMO')

%cal the C


