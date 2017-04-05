function [dz]=odefun(t,Z,Rot)

n=Rot.dim;
W = Rot.W;
A=Z(1:n);           %位移
B=Z((n+1):(n*2));   %速度

fd=zeros(n,1);
ufd=fd;    %不平衡力
bfd=fd;     %油腻压力

%计算 不平衡力
   ufx = Rot.RS.me*W^2*cos(W*t);
   %ufy = Rot.RS.me*W^2*sin(W*t);
   ufy = Rot.RS.me*W^2*cos(W*t)*W*t;
   
   %附值到ufd上;
   for bi = Rot.RS.Unban
       [unx,uny] = getXYnode(Rot,bi);
       ufd(unx) = ufx;
       ufd(uny) = ufy;
   end
   %计算 油膜压力
   for bi = Rot.RS.Springs
       [bnx,bny] = getXYnode(Rot,bi);
       [bfx,bfy] = getRotFxFyFun(Rot,A(bnx),A(bny),B(bnx),B(bny));
       %赋值到bfd上;
       bfd(bnx) = bfx;
       bfd(bny) = bfy;
   end
    fd = ufd+bfd;
    invM=inv(Rot.RS.mm);
    dz=zeros(2*n,1);
    dz(1:n) = B;
    dz((n+1):(n*2)) = (- invM*Rot.RS.cc*B-invM*Rot.RS.kk*A +fd);
end