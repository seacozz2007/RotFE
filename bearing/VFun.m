function [ V ] = VFun( x,y,a )
%�������Ĥ�����Ƽ��� ����V  VFUNG(x,y,a)
%   �˴���ʾ��ϸ˵��
V=(2+(y*cos(a)-x*sin(a))*GFun(x,y,a))/(1-x^2-y^2);
end

