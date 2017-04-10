    clear all;    
    options=odeset('RelTol',1e-3);   %acc
        init;
    [na nb]=size(M);
    iniM=zeros(1,nb*2);
    
    %iniM = iniM + 0.001;
    [t,xy]=ode45(@ffun,0:dt:dt*nt,iniM,options);

    for i=1:nb
        subplot(3,3,i+6); plot(t,xy(:,2*i-1),'-');
    end