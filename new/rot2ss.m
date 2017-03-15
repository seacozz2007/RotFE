function [a,b,c,d,bu,du]=rot2ss(varargin)
%[a,b,c,d,bu,du]=rot2ss(Rot)
% or
%[a,b,c,d,bu,du]=rot2ss(Rot,W)
% or
%[a,b,c,d,bu,du]=rot2ss(M,C,K,rc,fc)
% or
%[a,b,c,d,bu,du]=rot2ss('template'[,W])
%
%   
%Transform a rotor model (output of rotfe.m) to
% a state space. Optionally  reduce the order of the model
% according to the sub-field Reduct in Rot. If Rot.Reduct.flag==1
% the state transformation  Rot.T is being used to reduce the model
%  This transformation is computeddirectly by rotfe or by reduceX.m 
% 
%  The model described by this function is:
%
%   q'=a*q+b*u(t)+bu*(Fu_cos*cos(W*t)+Fu_sin*sin(W*t))
%   y=c*q
%   d,du are 0 and are provided for sole sake of compatability

% By  I. Bucher
% Date 15-3-1996   
% Rev. 1.0 
% 
% for Matlab 5.3  
% Rev  1.5  5-8-98
% Rev  2.01  2-8-99
input_form=3;
if isstruct(varargin{1}), % full FE etc. model in 
   switch nargin
   case 1,   Rot=varargin{1}; W=Rot.W;
   case 2, [Rot W]=deal(varargin{:});
   end
   input_form=1;
elseif isstr(varargin{1}) % rot2ss(File,...)
   Rot=rotfe( varargin{1} );
   W = Rot.W;
   if nargin>1, W=varargin{2}; end
   input_form=2;
end

switch input_form   
case {1,2}
   reduceQ=0;   
   if isfield(Rot,'Reduct'),
      if Rot.Reduct.flag==1, reduceQ=1;
      end,end
   
   if ~reduceQ % not reduced model >      
      n=size(Rot.M,1);
      if nargin>2, % dW(t)/dt-> alpha is provided
         a=[zeros(n) eye(n); -Rot.M\[Rot.K+C*Rot.Kst Rot.D+W*Rot.G]];
      else,
         a=[zeros(n) eye(n); -Rot.M\[Rot.K Rot.D+W*Rot.G]];
      end
      b=[zeros(n); inv(Rot.M)]; bu=b;
      c=[eye(n) zeros(n)];
      d=zeros(n);
      nr=length(Rot.FORCE_DOF); nf=length(Rot.RESP_DOF);
      if length(nr)>0, c=c(Rot.RESP_DOF,:); end,% only dof indicated in template
      if length(nf)>0,b=b(:,Rot.FORCE_DOF); end,% -"-
      if all([nf nr]>0), d=d(Rot.RESP_DOF,Rot.FORCE_DOF); end, % -"-
      du=zeros(size(c,1), size(bu,2));
      
   else
      n=size(Rot.M,1);
      a=[zeros(n) eye(n); -Rot.M\[Rot.K Rot.D+W*Rot.G]];
      bu=[zeros(n); inv(Rot.M)];
      b=bu*Rot.T.';
      d=zeros(n)*Rot.T.';
      
      c=[Rot.T 0*Rot.T];
      d=zeros(n)*Rot.T.';
      nr=length(Rot.FORCE_DOF); nf=length(Rot.RESP_DOF);   
      if nr>0, c=c(Rot.RESP_DOF,:); end, % only dof indicated in template
      if nf>0, b=b(:,Rot.FORCE_DOF); end, % -"-
      d=zeros( size(c,1), size(b,2)); % -"-
      du=zeros(size(c,1), size(bu,2));
   end % if ~reduceQ
   bu=full(bu); du=full(du);
   
case 3 % compatibility with former versions
   n=size(varargin{1},1);
   switch nargin
   case 5, [M,C,K,rc,fc]=deal(varargin{:});
   case 3, [M,C,K]=deal(varargin{:}); rc=1:n; fc=1:n;
   end
   
   
   a=[zeros(n) eye(n); -M\[K C]];
   b=[zeros(n); inv(M)];
   c=[eye(n) zeros(n)];
   d=zeros(n);
   b=b(:,fc); d=d(rc,fc); % take force cords & response
   c=c(rc,:) ;
end


a=full(a); b=full(b); c=full(c); d=full(d);

