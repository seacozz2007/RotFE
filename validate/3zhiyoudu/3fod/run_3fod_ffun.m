    clear all;    
    options=odeset('RelTol',1e-3);   %acc
    iniM=zeros(1,6);
    init_3fod;
    [t,xy]=ode45(@ffun_3fod,0:dt:dt*nt,iniM,options);
    subplot(3,3,7); plot(t,xy(:,1),'-');
    subplot(3,3,8); plot(t,xy(:,2),'-');
    subplot(3,3,9); plot(t,xy(:,3),'-');
