clear all;
%计算响应曲线
%尽量减少计算量


%global ufd bfd fd kk;
RSRot=rotfe('RotFE_model_1.m');
%建立初始模型

%建立响应计算参数
RotFE_Reduce_model_2;

RSRot.RS.T = 1;
RSRot.RS.Force=[];

RSRot.RS.Springs=[1 3];
RSRot.RS.Unban = [2];
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
%ode3doffun(t,Z,W,invM,kk,cc,me,S,C)
mm=RSRot.M;
kk=RSRot.K;
cc=RSRot.C;
invM=inv(M);
W=RSRot.W;
me=RSRot.RS.me;
initxy = zeros(12,1);
options=odeset('RelTol',1e-1);
[t,xy]=ode45(@ode3doffun,0:dt:nt*dt,initxy,options,W,invM,kk,cc,me,S,C);

subplot(3,1,1);
plot(t,xy(:,1));
subplot(3,1,2);
plot(t,xy(:,2));
subplot(3,1,3);
plot(t,xy(:,3));
