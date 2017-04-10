    clear all;    
    options=odeset('RelTol',1e-3);   %acc
    init;
    [na nb]=size(M);
    iniM=zeros(1,nb*2);

    %iniM = iniM + 0.001;
    [t,xy]=ode45(@ffun2,0:dt:dt*nt,iniM,options);
    for i=1:nb
        subplot(3,3,i+6); plot(t,xy(:,i),'-');
    end