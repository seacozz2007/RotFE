%计算响应曲线

global ufd bfd;
RSRot=rotfe('RotFE_model_1.m');
%建立初始模型

%建立响应计算参数
RotFE_model_2;

%计算响应曲线
[acc,vel,dsp]=fewilson(RSRot);


t = 0:RSRot.RS.dt:RSRot.RS.nt*RSRot.RS.dt;


plot(t,dsp(1,:));