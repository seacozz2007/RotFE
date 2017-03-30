function [Fx,Fy]=getRotFxFyFun(Rot,X,Y,X_dot,Y_dot)
    R = Rot.BEARING.R;
    C = Rot.BEARING.C;
    L = Rot.BEARING.L;
    u = Rot.BEARING.u;
    D = Rot.BEARING.D;
    x = X/C;
    y = Y/C;
    x_dot = X_dot/C;
    y_dot = Y_dot/C;
    S = u * Rot.W * (R^2/C^2)*(L^2/D^2)* R * L;
    [sFx,sFy] =  FxFyFun(x,y,x_dot,y_dot);
    Fx = sFx * S;
    Fy = sFy * S;
end