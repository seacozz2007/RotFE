clear all;
%������Ӧ����

%global ufd bfd fd kk;
RSRot=rotfe('RotFE_model_1.m');
%������ʼģ��

%������Ӧ�������
RotFE_model_2;

RSRot.RS.T = 1;
RSRot.RS.Force=[];

RSRot.RS.Springs=[];
RSRot.RS.Unban = [2];
%RSRot.RS.bcdof=[1 1 0 1 1 1 1 1 0 1 1 1];
%�ܲ���
RSRot.RS.nt = 1e3; 
%����
RSRot.RS.dt = 1e-1;
%ת��
RSRot.W = 1;
%�����ٻ�

%��������

%������Ӧ����

%������Ӧ����
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