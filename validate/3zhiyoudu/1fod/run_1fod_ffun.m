    clear all;    
    options=odeset('RelTol',1e-3);   %acc
    iniM=zeros(1,2);
    init_1fod;
    [t,xy]=ode45(@ffun_1fod,0:dt:dt*nt,iniM,options);
    subplot(2,1,2); plot(t,xy(:,1),'-');
