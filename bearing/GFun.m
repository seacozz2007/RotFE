function [ G ] = GFun( x,y,a )
%�������Ĥ�����Ƽ��� ����G  GFUNG(x,y,a)
%
    xy2=(1-x^2-y^2)^0.5;
    G=2/xy2*(pi/2+atan((y*cos(a)-x*sin(a))/xy2));
end