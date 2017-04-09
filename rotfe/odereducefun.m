function [dz]=odefun(t,Z,Rot)

n=Rot.dim;
W = Rot.W;
A=Z(1:n);           %λ��
B=Z((n+1):(n*2));   %�ٶ�

fd=zeros(n,1);
ufd=fd;    %��ƽ����
bfd=fd;     %����ѹ��

%���� ��ƽ����
   ufx = Rot.RS.me*W^2*cos(W*t);
   %ufy = Rot.RS.me*W^2*sin(W*t);
   ufy = Rot.RS.me*W^2*cos(W*t)*W*t;
   
   %��ֵ��ufd��;
   for bi = Rot.RS.Unban
       [unx,uny] = getXYnode(Rot,bi);
       ufd(unx) = ufx;
       ufd(uny) = ufy;
   end
   %���� ��Ĥѹ��
   for bi = Rot.RS.Springs
       [bnx,bny] = getXYnode(Rot,bi);
       [bfx,bfy] = getRotFxFyFun(Rot,A(bnx),A(bny),B(bnx),B(bny));
       %��ֵ��bfd��;
       bfd(bnx) = bfx;
       bfd(bny) = bfy;
   end
    fd = ufd+bfd;
    invM=inv(Rot.RS.mm);
    dz=zeros(2*n,1);
    dz(1:n) = B;
    dz((n+1):(n*2)) = (- invM*Rot.RS.cc*B-invM*Rot.RS.kk*A +fd);
end