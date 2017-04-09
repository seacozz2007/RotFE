function [nx,ny] = getXYnodeOde(Rot,n)
    nx = n;
    ny = n+Rot.dim/2;
end

