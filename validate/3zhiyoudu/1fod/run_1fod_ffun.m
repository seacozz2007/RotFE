    clear all;    
    options=odeset('RelTol',1e-3);   %acc
    iniM=zeros(1,2);
    init_1fod;
    [t,xy]=ode45(@ffun_1fod,0:dt:dt*nt,iniM,options);
    subplot(3,1,3); plot(t,xy(:,1),'-');
