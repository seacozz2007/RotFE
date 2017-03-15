Rot=rotfe('model_2node');
Rot.W=0;

KK=full(Rot.K);
MM=full(Rot.M);
N=2;
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
   NK(newrot,newrot)=eye(2);
   NM(newrot,newrot)=eye(2);
end

[NV,ND]=eig(NK,NM);
NW=sqrt(abs(diag(ND)));
NF=sort(NW)/2/pi

[VV,DD]=eig(KK,MM);
WW=sqrt(abs(diag(DD)));
FF=sort(WW)/2/pi

NE=eye(4);
NEK=mdiag(KK,NE);
NEM=mdiag(MM,NE);
[NEV,NED]=eig(NEK,NEM);
NEW=sqrt(abs(diag(NED)));
NEF=sort(NEW)/2/pi