clear all;
%计算响应曲线
%尽量减少计算量


%global ufd bfd fd kk;
RSRot=rotfe('RotFE_model_1.m');
%建立初始模型

%建立响应计算参数
RotFE_Reduce_model_2;

RSRot.RS.T = 1;
RSRot.RS.Force=[1];

RSRot.RS.Springs=[];
RSRot.RS.Unban = [];
%RSRot.RS.bcdof=[1 1 0 1 1 1 1 1 0 1 1 1];
%总步长
RSRot.RS.nt = 1e3; 
%步间
RSRot.RS.dt = 1e-1;

%无量纲化

R = RSRot.BEARING.R;
C = RSRot.BEARING.C;
L = RSRot.BEARING.L;
u = RSRot.BEARING.u;
D = RSRot.BEARING.D;

S = u * RSRot.W * (R^2/C^2)*(L^2/D^2)* R * L;
    

options=odeset();   %acc

nt=RSRot.RS.nt;
dt=RSRot.RS.dt;


% by:seaco 20170116
% m1=10;m2=10;m3=10;
% k1=1;k2=1;
% c1=0.1;c2=0.1;c3=0.1;
% 
% CM=[m1 0 0;0 m2 0;0 0 m3];
% CK=[k1 -k1 0;-k1 k1+k2 -k2;0 -k2 k2];
% CC=[c1 0 0;0 c2 0;0 0 c3];
% E0=zeros(3,3);
% RSRot.dim = 6;
% 
% RSRot.M=sparse([CM E0;E0 CM]);
% RSRot.K=sparse([CK E0;E0 CK]);
% RSRot.G=sparse([CC E0;E0 CC]);
% RSRot.C=sparse([CC E0;E0 CC]);


%ode3doffun(t,Z,W,invM,kk,cc,me,S,C)
mm=RSRot.M;
kk=RSRot.K;
cc=RSRot.C;

invM=inv(RSRot.M);
W=RSRot.W;
me=RSRot.RS.me;
initxy = zeros(12,1);
options=odeset('RelTol',1e-1);
[t,xy]=ode45(@ode3doffun,0:nt*dt,initxy,options,W,invM,kk,cc,me,S,C);

subplot(3,1,1);
plot(t,xy(:,1));
subplot(3,1,2);
plot(t,xy(:,2));
subplot(3,1,3);
plot(t,xy(:,3));
