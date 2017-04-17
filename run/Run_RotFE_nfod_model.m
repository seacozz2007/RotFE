% clear all;
% %������Ӧ����
% 
% %global ufd bfd fd kk;
% RSRot=rotfe('RotFE_model_1.m');
% %������ʼģ��
% 
% %������Ӧ�������
% RotFE_Reduce_model_2;
% 
% RSRot.RS.T = 1;
% RSRot.RS.Force=[];
% 
% RSRot.RS.Springs=[];
% RSRot.RS.Unban = [2];
% 
% %�ܲ���
% RSRot.RS.nt = 1e3; 
% %����
% RSRot.RS.dt = 1e-1;
% %ת��
% RSRot.W = 1;
% %�����ٻ�
% 
% %��������
% 
% 
% 
% 
% %������Ӧ����
% 

init_nfod_rotor;
[na n]=size(M);
kk=K;mm=M;cc=C;
q0=zeros(1,n);
dq0=zeros(1,n); 
bcdof=zeros(1,n);
fd=zeros(n,nt+1);

%���� fd
for it=1:(nt+1)                                              % loop for each time step
   %CAL THE F
    fd(1,it)=f1*cos(w*it*dt);
    %fd(:,it+1)=[f1*cos(w*(it+1)*dt)];
end

[acc,vel,dsp]=wilson(kk,cc,mm,fd,bcdof,nt,dt,q0,dq0);

t=0:dt:(dt*nt);
for i=1:n
    subplot(3,n,i); plot(t,dsp(i,:),'-');
end


%newmark

[acc,vel,dsp]=newmark(kk,cc,mm,fd,bcdof,nt,dt,q0,dq0);
t=0:dt:(dt*nt);
for i=1:n
    subplot(3,n,i+n); plot(t,dsp(i,:),'-');
end

   
    options=odeset('RelTol',1e-3);   %acc
    iniM=zeros(1,n*2);
    [t,xy]=ode45(@ffun_nfod_rotor,0:dt:dt*nt,iniM,options);
    for i=1:n
        subplot(3,n,i+n*2); plot(t,xy(:,i),'-');
    end