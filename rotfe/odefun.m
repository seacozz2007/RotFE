function [dz]=odefun(t,Z,Rot)

    n=Rot.dim;
    W = Rot.W;
    A=Z(1:n);           %位移
    B=Z((n+1):(n*2));   %速度

    fd=zeros(n,1);
    tfd=fd;     %激励力
    ufd=fd;    %不平衡力
    bfd=fd;     %油腻压力

    
    bcdof = Rot.RS.bcdof;
    

    %计算激励力 附值到tfd上;
    for bi = Rot.RS.Force
        tfx = Rot.RS.T*cos(W*t);
        tfy = Rot.RS.T*sin(W*t);
        [unx,uny] = getXYnode(Rot,bi);
        tfd(unx) = tfx;
        tfd(uny) = tfy;
    end


    %计算不平衡力 附值到ufd上;
    for bi = Rot.RS.Unban
        ufx = Rot.RS.me*W^2*cos(W*t);
        ufy = Rot.RS.me*W^2*cos(W*t)*W*t;
        [unx,uny] = getXYnode(Rot,bi);
        ufd(unx) = ufx;
        ufd(uny) = ufy;
    end
    %计算油膜压力 赋值到bfd上;
    for bi = Rot.RS.Springs
       [bnx,bny] = getXYnode(Rot,bi);
       [bfx,bfy] = getRotFxFyFun(Rot,A(bnx),A(bny),B(bnx),B(bny));
       
       bfd(bnx) = bfx;
       bfd(bny) = bfy;
    end
    fd = tfd + ufd+bfd;
    invM=inv(Rot.RS.mm);
    dz=zeros(2*n,1);
    dz(1:n) = B;
    dz((n+1):(n*2)) = (- invM*Rot.RS.cc*B-invM*Rot.RS.kk*A +fd);
    
    
     for i=1:n                  % assign zero to dsp, vel, acc of the dofs associated with bc
         if bcdof(i)==1
             dz(i)=0;
             dz(i+n)=0;
         end
     end
end