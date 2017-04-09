clear all;
%计算响应曲线

%global ufd bfd fd kk;
RSRot=rotfe('RotFE_model_1.m');
%建立初始模型

%建立响应计算参数
RotFE_model_2;

RSRot.RS.T = 1;
RSRot.RS.Force=[];

RSRot.RS.Springs=[];
RSRot.RS.Unban = [2];
%RSRot.RS.bcdof=[1 1 0 1 1 1 1 1 0 1 1 1];
%总步长
RSRot.RS.nt = 1e3; 
%步间
RSRot.RS.dt = 1e-1;
%转速
RSRot.W = 1;
%无量纲化

%参数矩阵

%计算响应曲线

%计算响应曲线
[acc,vel,dsp]=fewilson(RSRot);t2 = 0:RSRot.RS.dt:RSRot.RS.nt*RSRot.RS.dt;

subplot(2,2,1);
plot(t2,dsp(1,:));
subplot(2,2,2);
plot(t2,dsp(2,:));

% [t xy]=feode(RSRot);
% 
% subplot(2,2,3);
% plot(t,xy(:,1));
% subplot(2,2,4);
% plot(t,xy(:,2));