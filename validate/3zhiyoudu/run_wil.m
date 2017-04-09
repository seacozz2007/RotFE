clear all;
%振动方程组
%m1--k1--m2--k2--m3
%f1------f2------f3
%f1=sin(wt)  f2=0  f3=0
%c1------c2------c3

%方程:
%m1*x1..+c1*x1.+k1(x1-x2)=f1
%m2*x2..+c2*x2.+k1(x2-x1)+k2(x2-x3)=f2
%m3*x3..+c3*x3.+k2(x3-x2)=f3

m1=1;m2=1;m3=1;
k1=1;k2=1;
c1=0.1;c2=0.1;c3=0.1;
f1=1;
M=[m1 0 0;0 m2 0;0 0 m3];
K=[k1 -k1 0;-k1 k1+k2 -k2;0 -k2 k2];
C=[c1 0 0;0 c2 0;0 0 c3];

kk=K;mm=M;cc=C;
nt=1000;
dt=0.1;
q0=[0 0 0];
dq0=[0 0 0]; 
bcdof=[0 0 0];
%function [acc,vel,dsp]=wilson(kk,cc,mm,fd,bcdof,nt,dt,q0,dq0)

%--------------------------------------------------------------------------
%  Purpose:
%     The function subroutine TransResp4.m calculates transient response of
%     a structural system using Wilson   integration scheme.
%  Synopsis:
%     [acc, vel, dsp]=TransResp4(kk,cc,mm,fd,bcdof,nt,dt,q0,dq0)
%  Variable Description:
%     Input parameters
%       kk, cc, mm - stiffness, damping and mass matrices
%       fd - Input or forcing influence matrix
%       bcdof - Boundary condition dofs vector
%       nt - Number of time steps
%       dt - Time step size
%       q0, dq0 - Initial condition vectors
%     Output parameters
%       acc - Acceleration response
%       vel - Velocity response
%       dsp - Displacement response
%--------------------------------------------------------------------------
%  (1) initial condition
%--------------------------------------------------------------------------
[sdof,n2]=size(kk);

dsp=zeros(sdof,nt);                                         % displacement matrix
vel=zeros(sdof,nt);                                             % velocity matrix
acc=zeros(sdof,nt);                                          % acceleration matrix
 
dsp(:,1)=q0;                                               % initial displacement
vel(:,1)=dq0;                                                   % initial velocity
 
theta=1.4;                                          % select the parameters


%CAL THE F
fd(:,1)=[0 0 0];

%--------------------------------------------------------------------------
%  (2) Wilson   integration scheme for time integration
%--------------------------------------------------------------------------
acc(:,1)=inv(mm)*(fd(:,1)-kk*dsp(:,1)-cc*vel(:,1));
                                            % compute the initial acceleration (t=0)
ekk=kk+mm*(6/(theta*dt)^2)+cc*(3/(theta*dt));
                                           % compute the effective stiffness matrix
for i=1:sdof                  % assign zero to dsp, vel, acc of the dofs associated with bc
  if bcdof(i)==1
    dsp(i,1)=0;
    vel(i,1)=0;
    acc(i,1)=0;
  end
end

for it=1:nt                                              % loop for each time step
    
   %CAL THE F
    fd(:,it)=[sin(it*dt) 0 0];
    fd(:,it+1)=[sin((it+1)*dt) 0 0];
    
  cfm=dsp(:,it)*(6/(theta*dt)^2)+vel(:,it)*(6/(theta*dt))+2*acc(:,it);
  cfc=dsp(:,it)*(3/(theta*dt))+2*vel(:,it)+acc(:,it)*(theta*dt/2);
  efd=fd(:,it)+theta*(fd(:,it+1)-fd(:,it))+mm*cfm+cc*cfc;
                                            %  compute the effective force vector
  dtheta=inv(ekk)*efd;                         % find the displacement at time t+ dt

  acc(:,it+1)=(dtheta-dsp(:,it))*(6/(theta^3*dt^2))...
              -vel(:,it)*(6/(theta^2*dt))+acc(:,it)*(1-3/theta);
                                               % find the acceleration at time t+dt
  vel(:,it+1)=vel(:,it)+acc(:,it+1)*dt/2+acc(:,it)*dt/2;
                                                  % find the velocity at time t+dt
  dsp(:,it+1)=dsp(:,it)+vel(:,it)*dt...
              +(acc(:,it+1)+2*acc(:,it))*(dt^2/6);
                                              % find the displacement at time t+dt

  for i=1:sdof                % assign zero to acc, vel, dsp of the dofs associated with bc
    if bcdof(i)==1
      dsp(i,it+1)=0;
      vel(i,it+1)=0;
      acc(i,it+1)=0;
    end
  end
 
end
 
if cc(1,1)==0
  disp('The transient response results of undamping system')
else
  disp('The transient response results of damping system')
end

disp('The method is Wilson   integration')
%--------------------------------------------------------------------------
%     The end
%--------------------------------------------------------------------------
figure;
    t=0:dt:(dt*nt);
    subplot(3,1,1); plot(t,dsp(1,:),'-');
    subplot(3,1,2); plot(t,dsp(2,:),'-');
    subplot(3,1,3); plot(t,dsp(3,:),'-');
