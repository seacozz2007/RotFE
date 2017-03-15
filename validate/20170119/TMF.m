E=2e11;u=0.3;D1=0.02;p=7850;
A=pi*(D1/2)^2;
G=E/(2*(1+u));
Iz1=pi*D1^4/64;
l=0.1;
k1=3*E*Iz1/l^3;
pk=10/9;
%k1=G*A/pk/l;
M1=A*l*p;

M1=0.182755379347984;
M2=0.031929821979407;
k1=1.732495948670841e+07;
m=[M1/2 M1 M1 M1 M1 M1 M1 M1 M1 M1 M1/2];
k=[k1 k1 k1 k1 k1 k1 k1 k1 k1 k1];
    %TMM_TK方法求解固有频率
    %I为转动惯量，e为柔度
    N=length(m);
    M = diag(m);
    for i=1:(N-1);
       M(i,i+1) = M2; 
       M(i+1,i) = M2; 
    end
    K=zeros(N);
    for i=1:(N-1)
        K(i+1,i) = -k(i);
        K(i,i+1) = -k(i);
        K(i+1,i+1) =K(i+1,i+1) + k(i);
        K(i,i) = K(i,i) + k(i);
    end
    [V,D]=eig(K,M);

    %转化为列的形式，并求解出固有频率
    W=sqrt(abs(diag(D)));
    F = sort(W)/2/pi