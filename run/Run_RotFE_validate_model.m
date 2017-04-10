clear all;
%������Ӧ����
%�Ա�
%global ufd bfd fd kk;
%������ʼģ��
RSRot=rotfe('RotFE_model_1.m');

%������Ӧ�������
RotFE_validate_model_2;

%�ܲ���
RSRot.RS.nt = 1e4; 
%����
RSRot.RS.dt = 1e-1;
%ת��
RSRot.W = 1;
%�����ٻ�

%��������
RSRot.RS.kk = full(RSRot.K);
RSRot.RS.mm = full(RSRot.M);
RSRot.RS.gg = full(RSRot.G);
RSRot.invM = inv(RSRot.M);
RSRot.RS.cc = RSRot.RS.kk*RSRot.B+RSRot.W*RSRot.RS.gg;
%������Ӧ����
%[acc,vel,dsp]=fewilson(RSRot);t2 = 0:RSRot.RS.dt:RSRot.RS.nt*RSRot.RS.dt;

%plot(t,dsp(1,:));
%plot(t,fd(3,:));
%plot(t,dsp(1,:));

[t xy]=feode(RSRot);
subplot(2,2,1);
plot(t2,dsp(1,:));
subplot(2,2,2);
plot(t,xy(:,1));

subplot(2,2,3);
plot(t2,dsp(2,:));
subplot(2,2,4);
plot(t,xy(:,2));