% Non-linear dynamic analysis of a flexible rotor supported on porous oil journal bearings
% Nonlinear dynamic analysis of a flexible rotor supported by micropolar fluid film journal bearings
% script file which defines the model of a rotor, defines                                       
% (1) geometry                                                                                  
% (2) material                                                                                  
% (3) boundary conditions                                                                       
% (4) possible reduction of the model       
% (5) various flags which affect the run which follows     
% (6) Unbalance specification
% (7) Point force participation matrices
% (8) Point mass (linear and Angular, m & J)
%
% by:seaco 20170116

%阻尼系统
RSRot.B=25e-5;

%自由度
RSRot.RS.nDOF= RSRot.dim; 
%总步长
RSRot.RS.nt = 1000; 
%步间
RSRot.RS.dt = 0.1;

%激励力
RSRot.RS.Springs=[1 3];
%油膜压力
RSRot.RS.Springs=[1 3];
%不平衡力
RSRot.RS.Unban=[2];

%时间序列
RSRot.RS.t = 0:RSRot.RS.dt:(RSRot.RS.dt*RSRot.RS.nt);
%求解固有频率 转速为0
RSRot.FRes = eRotF(full(RSRot.K),full(RSRot.M));
%初始状态 位移
RSRot.RS.q0=zeros(RSRot.RS.nDOF,1);
%初始状态 速度
RSRot.RS.dq0=zeros(RSRot.RS.nDOF,1);
%边界条件  无约束
RSRot.RS.bcdof=zeros(RSRot.RS.nDOF,1);

%轴承参数
RSRot.BEARING.L = 0.01;
RSRot.BEARING.C = 90e-6;
RSRot.BEARING.D = 0.031;
RSRot.BEARING.R = RSRot.BEARING.D/2;
RSRot.BEARING.u = 0.04;
