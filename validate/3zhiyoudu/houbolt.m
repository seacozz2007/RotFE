function [acc,vel,dsp]=houbolt(kk,cc,mm,fd,bcdof,nt,dt,q0,dq0)
%--------------------------------------------------------------------------
%  Purpose:
%     The function subroutine TransResp3.m calculates transient response of
%     a structural system using Houbolt integration scheme.
%  Synopsis:
%     [acc,vel,dsp]=TransResp3(kk,cc,mm,fd,bcdof,nt,dt,q0,dq0)
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
 
dsp(:,1)=q0;                                         % initial displacement
vel(:,1)=dq0;                                            % initial velocity
%--------------------------------------------------------------------------
%  (2) calculate the initial integration using central difference scheme
%--------------------------------------------------------------------------
%---------------------------------------
%  (2.1) calculate the displacement at time -dt
%---------------------------------------
acc(:,1)=inv(mm)*(fd(:,1)-kk*dsp(:,1)-cc*vel(:,1));
                                           % compute the initial acceleration (t=0)
dsp0=dsp(:,1)-vel(:,1)*dt+0.5*acc(:,1)*dt^2;
                                    % compute the fictitious displacement at time -dt
 
ekk=mm/dt^2+cc/(2*dt);                       % compute the effective stiffness matrix
%----------------------------------------
%  (2.2) calculate the displacement at time t+dt
%----------------------------------------
it=1;
 
efd=fd(:,it)-(kk-2*mm/dt^2)*dsp(:,it)-(mm/dt^2-cc/(2*dt))*dsp0;
                                            %  compute the effective force vector
dsp(:,it+1)=inv(ekk)*efd;                                     % find the dsp at t+dt
                           % assign zero to acc, vel, dsp of the dofs associated with bc
for i=1:sdof
  if bcdof(i)==1
    dsp(i,1)=0;
    dsp(i,2)=0;
    vel(i,1)=0;
    acc(i,1)=0;
  end
end
%-----------------------------------------
%  (2.3) calculate the displacement at time t+2dt
%-----------------------------------------
for it=2:3                                       % loop for two steps after first step
  efd=fd(:,it)-(kk-2*mm/dt^2)*dsp(:,it)-(mm/dt^2-cc/(2*dt))*dsp(:,it-1);
                                            %  compute the effective force vector
  dsp(:,it+1)=inv(ekk)*efd;                                   % find the dsp at t+dt
 
  acc(:,it)=(dsp(:,it+1)-2*dsp(:,it)+dsp(:,it-1))/dt^2;                  % find the acc at t
  vel(:,it)=(dsp(:,it+1)-dsp(:,it-1))/(2*dt);                           % find the vel at t

  for i=1:sdof                  % assigh zero acc, vel, dsp of the dofs associated with bc
    if bcdof(i)==1
      dsp(i,it)=0;
      dsp(i,it+1)=0;
      vel(i,it)=0;
      acc(i,it)=0;
    end
  end
 
end
%--------------------------------------------------------------------------
%  (3) subsequent steps of the central difference scheme
%--------------------------------------------------------------------------
ekk=kk+cc*11/(6*dt)+mm*2/(dt*dt);
 
for it=3:nt                                 % loop for each time step after initial step
                                            %  compute the effective force vector
  cfm=dsp(:,it)*5/dt^2-dsp(:,it-1)*4/dt^2+dsp(:,it-2)*1/dt^2;
  cfc=dsp(:,it)*3/dt-dsp(:,it-1)*3/(2*dt)+dsp(:,it-2)*1/(3*dt);
 
  efd=fd(:,it+1)+mm*cfm+cc*cfc;
  
  dsp(:,it+1)=inv(ekk)*efd;                                   % find the dsp at t+dt
  acc(:,it)=(2*dsp(:,it+1)-5*dsp(:,it)+4*dsp(:,it-1)+dsp(:,it-2))/dt^2;
                                                          % find the acc at t+dt
  vel(:,it)=(11*dsp(:,it+1)-18*dsp(:,it)+9&dsp(:,it-1)-2*dsp(:,it-2))/(6*dt);
                                                          % find the vel at t+dt
  for i=1:sdof                % assign zero to acc, vel, dsp of the dofs associated with bc
    if bcdof(i)==1
      dsp(i,it)=0;
      dsp(i,it+1)=0;
      vel(i,it)=0;
      acc(i,it)=0;
    end
  end
 
end
 
if cc(1,1)==0
  disp('The transient response results of undamping system')
else
  disp('The transient response results of damping system')
end
disp('The method is Houbolt integration scheme')
%--------------------------------------------------------------------------
%     The end
% -------------------------------------------------------------------------
