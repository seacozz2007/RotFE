clear all;
%振动方程组
%m1--k1--m2--k2--m3
%f1------f2------f3
%f1=sin(wt)  f2=0  f3=0
%c1------c2------c3

%方程:
%m1*x1..+c1*x1.+k1(x1-x2)=f1
%m2*x2..+c2*x2.+k1(x2-x1)+k2(x2-x3)=f2
%m3*x3..+c3*x3.+k2(x3-x2)=f3

init;

kk=K;mm=M;cc=C;

[na nb]=size(kk);
q0=zeros(1,nb);
dq0=zeros(1,nb);
bcdof=zeros(1,nb);

fd=zeros(nb,nt+1);
%计算 fd
for it=1:(nt+1)                                              % loop for each time step
   %CAL THE F
    fd(1,it)=f1*cos(w*it*dt);
end

[acc,vel,dsp]=wilson(kk,cc,mm,fd,bcdof,nt,dt,q0,dq0);

    t=0:dt:(dt*nt);
for i=1:nb
        subplot(3,3,i); plot(t,dsp(i,:),'-');
    end