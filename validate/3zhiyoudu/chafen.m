function [acc,vel,dsp]=chafen(kk,cc,mm,fd,bcdof,nt,dt,q0,dq0)
%--------------------------------------------------------------------------
%  Purpose:
%     The function subroutine TransResp1.m calculates transient response of
%     a structural system using central difference scheme.
%  Synopsis:
%     [acc,vel,dsp]=TransResp1(kk,cc,mm,fd,bcdof,nt,dt,q0,dq0)
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

dsp(:,1)=q0;                                                % initial displacement
vel(:,1)=dq0;                                                   % initial velocity
%--------------------------------------------------------------------------
%  (2) central difference scheme for time integration
%--------------------------------------------------------------------------
acc(:,1)=inv(mm)*(fd(:,1)-kk*dsp(:,1)-cc*vel(:,1));
                                           % compute the initial acceleration (t=0)
dsp0=dsp(:,1)-vel(:,1)*dt+0.5*acc(:,1)*dt^2;
                                    % compute the fictitious displacement at time -dt
ekk=mm/dt^2+cc/(2*dt);
                                           % compute the effective stiffness matrix
%------------------------------------------
%  (2.1) first step of the central difference scheme
%------------------------------------------
efd=fd(:,1)-(kk-2*mm/dt^2)*dsp(:,1)-(mm/dt^2-cc/(2*dt))*dsp0;
                                            %  compute the effective force vector
dsp(:,1+1)=inv(ekk)*efd;
                                                         % find the dsp at 1+dt
for i=1:sdof                  % assign zero to acc, vel, dsp of the dofs associated with bc
  if bcdof(i)==1
    dsp(i,1)=0;
    dsp(i,2)=0;
    vel(i,1)=0;
    acc(i,1)=0;
  end
end
%-------------------------------------------------
%  (2.2) subsequent steps of the central difference scheme
%-------------------------------------------------
for it=2:nt                                   % loop for each time step after first step
  efd=fd(:,it)-(kk-2*mm/dt^2)*dsp(:,it)-(mm/dt^2-cc/(2*dt))*dsp(:,it-1);
                                            %  compute the effective force vector
  dsp(:,it+1)=inv(ekk)*efd;                                   % find the dsp at t+dt
  acc(:,it)=(dsp(:,it+1)-2*dsp(:,it)+dsp(:,it-1))/dt^2;                  % find the acc at t
  vel(:,it)=(dsp(:,it+1)-dsp(:,it-1))/(2*dt);                           % find the vel at t

  for i=1:sdof                % assign zero to acc, vel, dsp of the dofs associated with bc
    if bcdof(i)==1
      dsp(i,it)=0;
      dsp(i,it+1)=0;
      vel(i,it)=0;
      acc(i,it)=0;
    end
  end

end

  acc(:,it+1)=acc(:,it); vel(:,it+1)=vel(:,it);

if cc(1,1)==0
  disp('The transient response results of undamping system')
else
  disp('The transient response results of damping system')
end

disp('The method is central difference scheme 1')
%--------------------------------------------------------------------------
%     The end
%--------------------------------------------------------------------------
