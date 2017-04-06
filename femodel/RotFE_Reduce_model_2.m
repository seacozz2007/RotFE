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
MM = full(RSRot.M);
KK = full(RSRot.K);
GG = full(RSRot.G);

m1=MM(1,1);m2=MM(3,3);m3=MM(5,5);
c1=0.1;c2=0.1;c3=0.1;
f1=1;
M=[m1 0 0;0 m2 0;0 0 m3];
K=KK(1:2:5,1:2:5);
C=GG(1:2:5,1:2:5);

m1=1;m2=1;m3=1;
k1=1;k2=1;
c1=0.1;c2=0.1;c3=0.1;

%M=[m1 0 0;0 m2 0;0 0 m3];
K=[k1 -k1 0;-k1 k1+k2 -k2;0 -k2 k2];
C=[c1 0 0;0 c2 0;0 0 c3];

E0=zeros(3,3);
RSRot.dim = 6;

RSRot.M=sparse([M E0;E0 M]);
RSRot.K=sparse([K E0;E0 K]);
RSRot.G=sparse([C E0;E0 C]);
RSRot.C=sparse([C E0;E0 C]);
%���ɶ�
RSRot.RS.nDOF= RSRot.dim; 
%�ܲ���
RSRot.RS.nt = 1e3; 
%����
RSRot.RS.dt = 0.1;

%������
RSRot.RS.T = 1;
RSRot.RS.Force=[1];
%��Ĥѹ��
RSRot.RS.Springs=[];

%��ƽ����
RSRot.RS.me = 2e-5;
RSRot.RS.Unban=[];

%ʱ������
%RSRot.RS.t = 0:RSRot.RS.dt:(RSRot.RS.dt*RSRot.RS.nt);
%������Ƶ�� ת��Ϊ0

RSRot.FRes = eRotF(full(RSRot.K),full(RSRot.M));

%��ʼ״̬ λ��
RSRot.RS.q0=zeros(RSRot.RS.nDOF,1);
%��ʼ״̬ �ٶ�
RSRot.RS.dq0=zeros(RSRot.RS.nDOF,1);
%�߽�����  ��Լ��
RSRot.RS.bcdof=zeros(RSRot.RS.nDOF,1);
%RSRot.RS.bcdof=[1 1 0 0 1 1 1 1 0 0 1 1];

%��в���
RSRot.BEARING.L = 0.01;
RSRot.BEARING.C = 90e-6;
RSRot.BEARING.D = 0.031;
RSRot.BEARING.R = RSRot.BEARING.D/2;
RSRot.BEARING.u = 0.04;

%����ϵͳ
RSRot.B=0;
RSRot.W = 1;