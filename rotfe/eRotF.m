function [ FF ] = eRotF( K,M )
[V,D]=eig(K,M);
W=sqrt(abs(diag(D)));
FF=sort(W)/2/pi;
end

