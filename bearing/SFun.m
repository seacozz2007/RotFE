function [ S ] = SFun( x,y,a )
%短轴承油膜力近似计算 参数S  SFUNG(x,y,a)
%   
S=(x*cos(a)+y*sin(a))/(1-(x*cos(a)+y*sin(a))^2);
end

