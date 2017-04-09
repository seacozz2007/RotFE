%
% simply supported Euler beam
global L A Gs E Iy rho

L=1; E=2.15e11; rho=7800; d=50e-3; nu=0.3;
 r=d/2;

 ind=1:5;
  I=pi*d^4/64;   EI=E*I;
   A=pi*(r^2);
  lam=(ind*pi).^2/L^2*sqrt(EI/rho/A)/2/pi;

 
fprintf(' \n ');
%rotfe simply save
R=rotfe('simply');  R.W=0; R.MATERIALS(3)=0;
[v d]=roteig(R); w=(d)/2/pi; w=sort(w);
R=rotfe('simply');  R.W=0;  
[v d]=roteig(R); w1=(d)/2/pi; w1=sort(w1);


N=length(ind);
plot([w(1:2:2*N) w1(1:2:2*N) lam(1:N)'])


% use Genta page 87 eq. 2-40 to compute simply supported Timoshenko beam
% eigen frequencies

alpha=L*sqrt(A/I);  
xi=10/9; % in Genta's book
%xi=1.1282; % value in FE prog.
G=E/2/(1+nu); xis=xi*E/G;
 
Gs=G/xi;
Iy=I;
  

% solve determinant - d from MAPLE file timo_bm.ms
op=optimset('Display','off');
y=[]; for x=[100 600 2500 3000:2603:20000], y=[y fsolve('fs',x,op)]; end

% remove small ones
  y(find(y<20))=[];

  y=sort(y);   dy=diff(y); y(find(abs(dy)<1))=[];

y=y/2/pi;  % timoshenko beam frequencies

fprintf(' \n Euler Analytical ');
for q=1:length(lam),
   fprintf(' w%g=%8.5f ',q,lam(q));
 end
fprintf(' \n Euler FE ');
for q=1:length(lam),
   fprintf(' w%g=%8.5f ',q,w(2*q));
 end
fprintf(' \n Timoshenko FE ');

for q=1:length(lam),
   fprintf(' w%g=%8.5f ',q,w1(2*q));
 end
fprintf(' \n Timoshenko Analytical ');
 for q=1:length(y),
   fprintf(' w%g=%8.5f ',q,y(q));
 end
fprintf(' \n ');
fprintf(' \n Timoshenko difference FE-Analytical ');
 for q=1:length(y),
   fprintf(' w%g=%8.5f ',q,y(q)-w1(2*q));
 end
