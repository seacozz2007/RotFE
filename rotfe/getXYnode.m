function [nx,ny] = getXYnode(Rot,n)
    nx = n*2 - 1;
    ny = n*2-1+Rot.dim/2;
end

