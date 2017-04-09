function [Fx,Fy]=getSFxFyFun(S,x,y,x_dot,y_dot)
    [sFx,sFy] =  FxFyFun(x,y,x_dot,y_dot);
    Fx = sFx * S;
    Fy = sFy * S;
end