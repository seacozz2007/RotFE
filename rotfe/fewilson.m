function [acc,vel,dsp]=fewilson(Rot)
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
kk=full(Rot.K);
mm=full(Rot.M);
cc=Rot*kk*Rot.B+Rot.W*full(Rot.G);
[sdof,n2]=size(kk);
nt=Rot.RS.nt;
dsp=zeros(sdof,nt);                                         % displacement matrix
vel=zeros(sdof,nt);                                             % velocity matrix
acc=zeros(sdof,nt);                                          % acceleration matrix
 
dsp(:,1)=Rot.RS.q0;                                               % initial displacement
vel(:,1)=Rot.RS.dq0;                                                   % initial velocity
W = RS.W;
bcdof = Rot.RS.bcdof;
theta=1.4;                                          % select the parameters
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
fd=zeros(RS.dim,nt);
for it=1:nt                                              % loop for each time step
   %计算 不平衡力
   
   %x方向
   fx = Rot.RS.me*W^2*cos(W*it*dt);
   %fy = Rot.RS.me*W^2*sin(W*it*dt);
   fy = Rot.RS.me*W^2*cos(W*it*dt)*W*it*dt;
   
   %计算 油膜压力
   
   fd(Rot.RS.Unban*2-1,it+1) = fx;
   fd(Rot.RS.Unban*2-1+Rot.dim/2) =fy;
   %附值到fd上;
   
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
