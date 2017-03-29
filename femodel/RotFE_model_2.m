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

%����ϵͳ
RSRot.B=25e-5;

%���ɶ�
RSRot.RS.nDOF= RSRot.dim; 
%�ܲ���
RSRot.RS.nt = 1000; 
%����
RSRot.RS.dt = 0.1;

%������
RSRot.RS.Springs=[1 3];
%��Ĥѹ��
RSRot.RS.Springs=[1 3];
%��ƽ����
RSRot.RS.Unban=[2];

%ʱ������
RSRot.RS.t = 0:RSRot.RS.dt:(RSRot.RS.dt*RSRot.RS.nt);
%������Ƶ�� ת��Ϊ0
RSRot.FRes = eRotF(full(RSRot.K),full(RSRot.M));
%��ʼ״̬ λ��
RSRot.RS.q0=zeros(RSRot.RS.nDOF,1);
%��ʼ״̬ �ٶ�
RSRot.RS.dq0=zeros(RSRot.RS.nDOF,1);
%�߽�����  ��Լ��
RSRot.RS.bcdof=zeros(RSRot.RS.nDOF,1);

%��в���
RSRot.BEARING.L = 0.01;
RSRot.BEARING.C = 90e-6;
RSRot.BEARING.D = 0.031;
RSRot.BEARING.R = RSRot.BEARING.D/2;
RSRot.BEARING.u = 0.04;
