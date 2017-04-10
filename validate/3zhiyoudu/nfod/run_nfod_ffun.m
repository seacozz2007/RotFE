    clear all;    
    options=odeset('RelTol',1e-3);   %acc
    init_nfod;
    iniM=zeros(1,n*2);
    [t,xy]=ode45(@ffun_nfod,0:dt:dt*nt,iniM,options);
    for i=1:n
        subplot(3,n,i+n*2); plot(t,xy(:,i),'-');
    end
