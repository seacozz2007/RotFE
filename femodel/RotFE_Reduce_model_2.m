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

M=MM(1:2:12,1:2:12);
%��������һ����ȷ;
%M=diag(diag(M),0);

K=KK(1:2:12,1:2:12);
G=GG(1:2:12,1:2:12);


RSRot.dim = 6;

RSRot.M=sparse(M);
RSRot.K=sparse(K);
RSRot.G=sparse(G);

%���ɶ�
RSRot.RS.nDOF= RSRot.dim; 
%�ܲ���
RSRot.RS.nt = 1e3; 
%����
RSRot.RS.dt = 0.1;

%������
RSRot.RS.T = 1;
RSRot.RS.Force=[];
%��Ĥѹ��
RSRot.RS.Springs=[1 3];

%��ƽ����
RSRot.RS.me = 2e-5;
RSRot.RS.Unban=[2];

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
RSRot.B=25e-5;
RSRot.W = 20*pi;