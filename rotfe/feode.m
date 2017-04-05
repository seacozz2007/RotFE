function [t,xy] = feode(Rot)
    options=odeset('RelTol',1e-2);   %acc
    nt=Rot.RS.nt;
    dt=Rot.RS.dt;
    [t,xy]=ode45(@odefun,[0:dt:nt*dt],[Rot.RS.q0 Rot.RS.dq0],options,Rot);
    %figure
    %subplot(2,1,1);
    %plot(t(:),xy(:,1),'-'); 
    %subplot(2,1,2);
    %plot(t(:),xy(:,3),'-');

end