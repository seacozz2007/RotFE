function [dz]=ode3doffun(t,Z,W,invM,kk,cc,me,S,C)
    t=t
    n=6;
    %1:x1 2:x2 3:x3 4:y1 5:y2 6:y3

    A=Z(1:n);           %位移
    B=Z((n+1):(n*2));   %速度

    fd=zeros(n,1);


    %ufx = me*W^2*cos(W*t);
    %ufy = me*W^2*sin(W*t);
    %ufy = me*W^2*cos(W*t)*W*t;
   
    fd(2) = me*W^2*cos(W*t);
    fd(5) = me*W^2*cos(W*t)*W*t;

    %计算 油膜压力
%     [bfx,bfy] = getSFxFyFun(S,A(1)/C,A(4)/C,B(1)/C,B(4)/C);
%     fd(1) = bfx;
%     fd(4) = bfy;
%     [bfx,bfy] = getSFxFyFun(S,A(3)/C,A(6)/C,B(3)/C,B(6)/C);
%     fd(3) = bfx;
%     fd(6) = bfy;

    dz=zeros(2*n,1);
    dz(1:n) = B;
    dz((n+1):(n*2)) = (- invM*cc*B-invM*kk*A + fd);
end