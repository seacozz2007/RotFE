function [ A ] = AFun( x,y,x_dot,y_dot )
%AFUN �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
if((x-2*y_dot) == 0)
    A = 0;
    return;
end
tmp_A = (y+2*x_dot)/(x-2*y_dot);
A=atan(tmp_A)-pi/2*sign(tmp_A)-pi/2*sign(y+2*x_dot);

end

