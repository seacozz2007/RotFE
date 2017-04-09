    clear all;    
    options=odeset('RelTol',1e-3);   %acc
    iniM=zeros(1,6);
    
    %iniM = iniM + 0.001;
    [t,xy]=ode45(@ffun,0:0.1:100,iniM,options);
    subplot(3,3,4); plot(t,xy(:,1),'-');
    subplot(3,3,5); plot(t,xy(:,3),'-');
    subplot(3,3,6); plot(t,xy(:,5),'-');