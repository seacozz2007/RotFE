    clear all;    
    options=odeset('RelTol',1e-3);   %acc
    iniM=zeros(1,6);

    %iniM = iniM + 0.001;
    [t,xy]=ode45(@ffun2,0:0.1:100,iniM,options);
    subplot(3,3,7); plot(t,xy(:,1),'-');
    subplot(3,3,8); plot(t,xy(:,2),'-');
    subplot(3,3,9); plot(t,xy(:,3),'-');