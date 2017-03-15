R21=rotfe('asym21');


[r s]=roteig1(R21);

q=12;
[Wx Wy]=rotrev(R21,r(:,q));
h=plot3(R21.NODES/5,Wx(1:2:end,:),Wy(1:2:end,:),'-.'); % plot linear displacements only

set(h(1),'linewidth',4)