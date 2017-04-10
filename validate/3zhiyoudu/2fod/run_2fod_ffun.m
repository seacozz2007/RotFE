    clear all;    
    options=odeset('RelTol',1e-3);   %acc
    iniM=zeros(1,4);
    init_2fod;
    [t,xy]=ode45(@ffun_2fod,0:dt:dt*nt,iniM,options);
    subplot(3,2,5); plot(t,xy(:,1),'-');
    subplot(3,2,6); plot(t,xy(:,2),'-');
