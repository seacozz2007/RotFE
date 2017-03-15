Rot=rotfe('m1');

Rot.W=0;
[phi lambda]=roteig(Rot);
Rot.eigenvector=phi; Rot.eigenvalue=abs(lambda); 
F=abs(lambda)/2/pi
