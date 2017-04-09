function [t,xy] = feode(Rot)
   
    options=odeset();   %acc
     %options=odeset('RelTol',1e-3);   %acc
    
    nt=Rot.RS.nt;
    dt=Rot.RS.dt;
    [t,xy]=ode45(@odefun,0:nt*dt,[Rot.RS.q0 Rot.RS.dq0],options,Rot);
end