%
% simply supported Euler beam
global L A Gs E Iy rho nn
% use Genta page 87 eq. 2-40 to compute simply supported Timoshenko beam
% eigen frequencies
dw=[];
%
L=1; E=2.15e11; rho=7800; d=50e-3; nu=0.3;
r=d/2;

ind=1:5;
I=pi*d^4/64;   EI=E*I;
A=pi*(r^2);
lam=(ind*pi).^2/L^2*sqrt(EI/rho/A)/2/pi;

alpha=L*sqrt(A/I);  
xi=10/9; % in Genta's book
xi=1.1282; % value in FE prog. (more accurate)
G=E/2/(1+nu); xis=xi*E/G;

Gs=G/xi;
Iy=I;


% solve determinant - d from MAPLE file timo_bm.ms
op=optimset('Display','off','TolX',1e-4,'TolFun',1e-4);
y=[]; for x=[100 600 2500 3000:2603:20000], y=[y fsolve('fs',x,op)]; end

% remove small ones
y(find(y<20))=[];

y=sort(y);   dy=diff(y); y(find(abs(dy)<1))=[];

y=y/2/pi;  % timoshenko beam frequencies


for nn=[1 2 5 10 20 40 ]
   
   fprintf(' \n ');
   %rotfe simply save
   R=rotfe('simplyx');  R.W=0; R.MATERIALS(3)=0;   R=rotfe(R); 
   [v d]=roteig(R); w=abs(d)/2/pi; w=sort(w);
   R=rotfe('simplyx');  R.W=0;  
   [v d]=roteig(R); w1=abs(d)/2/pi; w1=sort(w1);
   
   
   N=length(ind);
   %plot([w(1:2:2*N) w1(1:2:2*N) lam(1:N)'])
   
   
   fprintf(' \n Euler Analytical ');
   for q=1:length(lam),
      fprintf(' w%g=%8.5f ',q,lam(q));
   end
   fprintf(' \n Euler FE ');
   for q=1:length(w)/2,
      fprintf(' w%g=%8.5f ',q,w(2*q));
   end
   fprintf(' \n Timoshenko FE ');
   
   for q=1:length(w1)/2,
      fprintf(' w%g=%8.5f ',q,w1(2*q));
   end
   fprintf(' \n Timoshenko Analytical ');
   for q=1:length(y),
      fprintf(' w%g=%8.5f ',q,y(q));
   end
   fprintf(' \n ');
   fprintf(' \n Timoshenko difference FE-Analytical ');
   for q=1:min(length(y),length(w1)/2),
      fprintf(' w%g=%8.5f ',q,y(q)-w1(2*q));
   end
   n1=min(length(y),length(w1)/2);
   dww=[y(1:n1)'-w1(1:2:n1*2)]; if size(dww,1)<5, dww(5)=0; end
   dw=[dw  dww];
end
loglog( [1 2 5 10 20 40  ]'*(dw(:,1)*0+1)', abs(dw)','-.')
legend('w1','w2','w3','w4','w5')
title(' analytic - FE, 1st 5 natural frequencies')
xlabel(' no. of elements')
ylabel('error')