function T=get_view()
%
% same as T=view; but without Matlab's bug :-) !!!

 T=view;   T(3,1:3) = -T(3,1:3);
