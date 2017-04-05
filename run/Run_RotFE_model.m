clear all;
%������Ӧ����

%global ufd bfd fd kk;
RSRot=rotfe('RotFE_model_1.m');
%������ʼģ��

%������Ӧ�������
RotFE_model_2;

RSRot.RS.T = 1;
RSRot.RS.Force=[2];

RSRot.RS.Springs=[];
RSRot.RS.Unban = [];
%RSRot.RS.bcdof=[1 1 0 1 1 1 1 1 0 1 1 1];
%�ܲ���
RSRot.RS.nt = 1e2; 
%����
RSRot.RS.dt = 1e-1;
%ת��
RSRot.W = 1;
%�����ٻ�

%��������
RSRot.RS.kk = full(RSRot.K);
RSRot.RS.mm = full(RSRot.M);
RSRot.RS.gg = full(RSRot.G);
RSRot.RS.cc=RSRot.RS.kk*RSRot.B+RSRot.W*RSRot.RS.gg;
%������Ӧ����
[acc,vel,dsp]=fewilson(RSRot);t2 = 0:RSRot.RS.dt:RSRot.RS.nt*RSRot.RS.dt;

%plot(t,dsp(1,:));
%plot(t,fd(3,:));
%plot(t,dsp(1,:));

[t xy]=feode(RSRot);
subplot(2,1,1);
plot(t2,dsp(3,:));
subplot(2,1,2);
plot(t,xy(:,3));