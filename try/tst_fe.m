

 d=3e-3; L=.16;
 E=210e9; rho=7800;
 F=2*9.81;
 X=[];
for F=[0:.25:4]*98.1

 n=59;
 ne=2*(n+1);

 M=zeros(ne); K=zeros(ne);

 for q=1:n,

  [Me,Ke]=mkbeam(d,L/n,E,rho,-F);

 ind=(q-1)*2+1:(q-1)*2+4;
 M(ind,ind)= M(ind,ind)+Me;
 K(ind,ind)= K(ind,ind)+Ke;

end

 ind=[1 2 ne-1 ne];
K(ind,:)=[]; K(:,ind)=[];
M(ind,:)=[]; M(:,ind)=[];

 [vv dd]=eig(K,M);

 plot([0;vv(1:2:end,1);0])
title(sprintf('f=%g',sqrt(dd(1))/2/pi));
drawnow
 hold on
 X=[X; sqrt(dd(1))/2/pi];
end




