function [ V ] = VFun( x,y,a )
%短轴承油膜力近似计算 参数V  VFUNG(x,y,a)
%   此处显示详细说明
V=(2+(y*cos(a)-x*sin(a))*GFun(x,y,a))/(1-x^2-y^2);
end

