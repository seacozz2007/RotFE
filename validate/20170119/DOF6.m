Rot=rotfe('m2');
Rot.W=0;

KK=full(Rot.K);
MM=full(Rot.M);
N=11;
oldDOF=4;
newDOF=6;

NK=zeros(newDOF*N);
NM=zeros(newDOF*N);
for i=1:N
   oldse=(i-1)*oldDOF+1:i*oldDOF;
   newse=(i-1)*newDOF+1:(i-1)*newDOF+oldDOF;
   
   newrot=((i-1)*newDOF+oldDOF+1):i*newDOF;
   NK(newse,newse)=KK(oldse,oldse);
   NM(newse,newse)=MM(oldse,oldse);
   
   %rotNK
   NK(newrot,newrot)=[1,1;1,1];
   NM(newrot,newrot)=[1,1;1,1];
end

[NV,ND]=eig(NK,NM);
NW=sqrt(abs(diag(ND)));
NF=sort(NW)/2/pi

[VV,DD]=eig(KK,MM);
WW=sqrt(abs(diag(DD)));
FF=sort(WW)/2/pi

NE=zeros(10);
KK=mdiag(KK,NE);
MM=mdiag(MM,NE);
[VV,DD]=eig(KK,MM);
WW=sqrt(abs(diag(DD)));
FF=sort(WW)/2/pi