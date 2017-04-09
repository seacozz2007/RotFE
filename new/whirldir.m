function [ws,wd]=whirldir(v)
% ws=whirldir(v)
%
% find direction of whirl
%  ws=1,	forward
%  ws=-1,	backward
%  ws=0,		standing wave or mixed mode (forward and backward)
%
% by I. Bucher  28.4.1995
%
% v is a complex eigenvector
%
% v = [x ; y]  (or [y ; z])

[n m1]=size(v);
n1=fix(n/2);
if 2*n1~=n,	
   errordlg('ERROR in whirldir',' length(v) should be even','modal'); 
   return, 
end


for q=1:m1,
   
   a=v(1:n1,q); b=v(n1+1:n,q)/i;	        % separate (see whirldir.ms, MAPLE file)
   s=(real(a).*real(b)+imag(a).*imag(b)); % that's the inequality
   ns=norm(s)/length(s);
   if all(ns==0), ns=1; end
   s( find( abs(s/ns)<1e-3) )=[]; % remove small ones
%   s(find(abs(s)<(norm(s)/length(s)*1e-4)))=[];
%   TOL=3e-4*norm(s)/sqrt([norm(a)^2+norm(b)^2]);
%   s=s.*(abs(s)>TOL);
   s=(s>0)-(s<0); 
   
   na=norm(a)+eps;  nb=norm(b)+eps;  SWR=na/nb+nb/na;
   if (all(s==0)) | SWR>300, ws(q)=0; wd=' Standing wave';
   elseif all(s<0),  ws(q)=1; wd=' Forward'; 				% forward
   elseif all(s>0), ws(q)=-1;wd=' Backward'; 				% backward
   else,ws(q)=0; wd=' Mixed Forward/Backward';end;					% undefined
   
end
 
