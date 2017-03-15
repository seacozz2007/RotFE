%build the fe model
Rot=rotfe('ROT')

%cal the C
b=25e-5;
w=20;

Rot.C=b*Rot.K+w*Rot.G;
