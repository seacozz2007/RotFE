function [acc,vel,dsp]=newmark(kk,cc,mm,fd,bcdof,nt,dt,q0,dq0)
%--------------------------------------------------------------------------
%  Purpose:
%     The function subroutine TransResp5.m calculates transient response of
%     a structural system using Newmark integration scheme.
%  Synopsis:
%     [acc,vel,dsp]=TransResp5(kk,cc,mm,fd,bcdof,nt,dt,q0,dq0)
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
vel(:,1)=dq0;                                                  % initial velocity

alpha=0.5; beta=0.5;                                        % select the parameters
%--------------------------------------------------------------------------
%  (2) Newmark integration scheme for time integration
%--------------------------------------------------------------------------
acc(:,1)=inv(mm)*(fd(:,1)-kk*dsp(:,1)-cc*vel(:,1));
                                           % compute the initial acceleration (t=0)
ekk=kk+mm/(alpha*dt^2)+cc*beta/(alpha*dt);
                                          % compute the effective stiffness matrix
                           % assign zero to dsp, vel, acc of the dofs associated with bc
for i=1:sdof
  if bcdof(i)==1
    dsp(i,1)=0;
    vel(i,1)=0;
    acc(i,1)=0;
  end
end

for it=1:nt                                              % loop for each time step
  cfm=dsp(:,it)/(alpha*dt^2)+vel(:,it)/(alpha*dt)+acc(:,it)*(0.5/alpha-1);
  cfc=dsp(:,it)*beta/(alpha*dt)+vel(:,it)*(beta/alpha-1)...
     +acc(:,it)*(0.5*beta/alpha-1)*dt;
  efd=fd(:,it)+mm*cfm+cc*cfc;                  %  compute the effective force vector

  dsp(:,it+1)=inv(ekk)*efd;                        % find the displacement at time t+dt
  acc(:,it+1)=(dsp(:,it+1)-dsp(:,it))/(alpha*dt^2)-vel(:,it)/(alpha*dt)...
              -acc(:,it)*(0.5/alpha-1);              % find the acceleration at time t+dt
  vel(:,it+1)=vel(:,it)+acc(:,it)*(1-beta)*dt+acc(:,it+1)*beta*dt;
                                                  % find the velocity at time t+dt
                           % assign zero to acc, vel, dsp of the dofs associated with bc
  for i=1:sdof
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
 
disp('The method is Newmark integration')
%--------------------------------------------------------------------------
%     The end
%--------------------------------------------------------------------------
