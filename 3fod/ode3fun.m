function [dz]=ode3fun(t,Z)
    
    W = 1;
    m1 = 1;m2 = 1;m3 = 1;
    k1 = 1;k2 = 1;
    c1 = 1;c2 = 1;c3 = 1;
    n = 3;
    
    dz=zeros(12,1);
    %x1: Z(1) x1'=Z(2)
    dz(1) = Z(2);
    dz(2) = c1  

%     %���㼤���� ��ֵ��tfd��;
%     for bi = Rot.RS.Force
%         tfx = Rot.RS.T*cos(W*t);
%         tfy = Rot.RS.T*sin(W*t);
%         [unx,uny] = getXYnodeOde(Rot,bi);
%         tfd(unx) = tfx;
%         tfd(uny) = tfy;
%     end


    %���㲻ƽ���� ��ֵ��ufd��;

        ufx = me*W^2*cos(W*t);
        %ufy = me*W^2*cos(W*t)*W*t;
        
        fd(2) = me*W^2*cos(W*t);

    
    %������Ĥѹ�� ��ֵ��bfd��;
    for bi = Rot.RS.Springs
       [bnx,bny] = getXYnodeOde(Rot,bi);
       [bfx,bfy] = getRotFxFyFun(Rot,A(bnx),A(bny),B(bnx),B(bny));
       
       bfd(bnx) = bfx;
       bfd(bny) = bfy;
    end
    fd = tfd + ufd + bfd;
    
    invM=inv(Rot.M);
    
    dz=zeros(2*n,1);
    dz(1:n) = B;
    dz((n+1):(n*2)) = (- invM*Rot.C*B-invM*Rot.K*A +invM*fd);
        
    for i=1:n                  % assign zero to dsp, vel, acc of the dofs associated with bc
         if bcdof(i)==1
             dz(i)=0;
             dz(i+n)=0;
         end
     end
end