E=2e11;pa=0.3;D1=0.02;D2=0.03;p=7850;
Iz1=pi*D1^4/64;Iz2=pi*D2^4/64;
l=0.1;
k1=3*E*Iz1/l^3;k2=3*E*Iz2/l^3;
M1=pi*(D1/2)^2*l*p;
M2=pi*(D2/2)^2*l*p;

m=[M1/2 M1 M1 M1 M1 (M1/2+M2/2) M2 M2 M2 M2];

    %TMM_TK方法求解固有频率
    %I为转动惯量，e为柔度
    N=length(m);
    M = zeros(N);
    for i=1:N
       M(i,i) = m(i); 
    end
    K=zeros(N);
    for i=1:(N-1)
        K(i+1,i) = -k(i);
        K(i,i+1) = -k(i);
        K(i+1,i+1) =K(i+1,i+1) + k(i);
        K(i,i) = K(i,i) + k(i);
    end
    [V,D]=eig(K,M);
    W=zeros(N,1);
    %转化为列的形式，并求解出固有频率
    for i=1:N
    W(i)=sqrt(D(i,i));
    end
    W;
    F = W/2/pi;