%������Ӧ����

global ufd bfd;
RSRot=rotfe('RotFE_model_1.m');
%������ʼģ��

%������Ӧ�������
RotFE_model_2;

%������Ӧ����
[acc,vel,dsp]=fewilson(RSRot);


t = 0:RSRot.RS.dt:RSRot.RS.nt*RSRot.RS.dt;


plot(t,dsp(1,:));