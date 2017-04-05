function [acc,vel,dsp]=fewilson(Rot)
%global ufd bfd fd kk;
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
cc=kk*Rot.B+Rot.W*full(Rot.G);
[sdof,n2]=size(kk);
nt=Rot.RS.nt;
dt=Rot.RS.dt;

dsp=zeros(sdof,nt);                                         % displacement matrix
vel=zeros(sdof,nt);                                             % velocity matrix
acc=zeros(sdof,nt);                                          % acceleration matrix

fd=zeros(Rot.dim,nt+1);
tfd=fd;     %激励力
ufd=fd;     %不平衡力
bfd=fd;     %油腻压力

dsp(:,1)=Rot.RS.q0;                                               % initial displacement
vel(:,1)=Rot.RS.dq0;                                                   % initial velocity
W = Rot.W;
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

for it=1:nt                                              % loop for each time step
    
   %计算激励力 附值到tfd上;
    for bi = Rot.RS.Force
        tfx = Rot.RS.T*cos(W*it*dt);
        tfy = Rot.RS.T*sin(W*it*dt);
        [unx,uny] = getXYnode(Rot,bi);
        tfd(unx,it+1) = tfx;
        tfd(uny,it+1) = tfy;
    end
    
   
    %计算不平衡力 附值到ufd上;
    for bi = Rot.RS.Unban
        ufx = Rot.RS.me*W^2*cos(W*it*dt);
        ufy = Rot.RS.me*W^2*cos(W*it*dt)*W*it*dt;
        [unx,uny] = getXYnode(Rot,bi);
        ufd(unx,it+1) = ufx;
        ufd(uny,it+1) = ufy;
    end

    %计算 油膜压力
    for bi = Rot.RS.Springs
        [bnx,bny] = getXYnode(Rot,bi);
        [bfx,bfy] = getRotFxFyFun(Rot,dsp(bnx,it),dsp(bny,it),vel(bnx,it),vel(bny,it));
        %赋值到bfd上;
        bfd(bnx,it+1) = bfx;
        bfd(bny,it+1) = bfy;
    end
    fd(:,it+1) = tfd(:,it+1) + ufd(:,it+1)+bfd(:,it+1);

  %display([' x:' num2str(dsp(bnx,it)) '   y:' num2str(dsp(bny,it))])
  %display(['cx:' num2str(dsp(bnx,it)/Rot.BEARING.C) '  cy:' num2str(dsp(bny,it)/Rot.BEARING.C)])
  %display(['ufx:' num2str(ufx) ' ufy:' num2str(ufy)])
  %display(['bfx:' num2str(bfx) ' bfy:' num2str(bfy)])
  %display([' fx:' num2str(bfx) '  fy:' num2str(bfy)])
  
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

