n=2;
m1=4;m2=7;m3=3;
k1=16;k2=3;k3=3;k4=3;
c1=0.3;c2=0.3;c3=0.2;
f1=1;
w=5;

M=[m1 0 0;0 m2 0;0 0 m3];
K=[k1+k2 -k2 0;-k2 k2+k3 -k3;0 -k3 k3+k4];
C=[c1 0 0;0 c2 0;0 0 c3];

nt=10000;
dt=0.01;