function [phio,lambda,fb]=rotmodescompute(file,varargin) 
% [phi,lambda,fb]=rotmodes(file,W,TYPE,Nmodes,TOL) 
%   
% compute (for a rotor): 
% (1) right hand modes
% (2) eigenvalues
% (3) sense of rotation (forward/backward or mixed)
%
% INPUTS:
%	FILE:  input file, can be template.m (Script) or template.mtx *.mat
%	W:	speed of rotation [Rad/s] (default is W=0)
%	TYPE:	type of eigenproblem to be solved
%		TYPE=0,  undamped (default)
%		TYPE=1,  external (non-rotating) damping
%		TYPE=2,  general (rotating hysteretic + non-rotating viscous) damping
% EXAMPLE:
%	 [phi,lambda,fb]=rotmodes('template');  % runs rotfe('template'), 
%						% assumes W=0

% By  I. Bucher  8-4-1996
% Rev. 1.01 for RRA   

% process inputs, find any command strings
% currently supported:
%	save,  plot

% define other inputs if not provided
 Vars=varargin;  
  if nargin<2, W=0;       else, W=Vars{1};      end, if isstr(W), W=0; end
  if nargin<3, TYPE=0;    else, TYPE=Vars{2};   end, if isstr(TYPE), TYPE=0; end
  if nargin<4, Nmodes=10; else, Nmodes=Vars{3}; end, if isstr(Nmodes), Nmodes=10; end
  if nargin<5, TOL=1e-3;  else, TOL=Vars{4};    end, if isstr(TOL), TOL=1e-3; end
  cmp=['save' ; 'plot']; OP=[0 ; 0];
  for q=1:nargin-1,
    for p=1:size(cmp,1),
       OP(p)=OP(p) | strcmp(Vars{q}, cmp(p,:));
    end
  end         	
Save=OP(1); Plot=OP(2);
Ws=W;

% process file *.m or *.mat file
if isstr(file),
 FAIL=['error([''rotmodes:  file >> '' FILE '' << not found in MATLAB path''])'];
 [FILE,ext]=nameext(file);	% get name without extention
 if strcmp(ext,'mat'),
	cmd=['load ' FILE]; runQ=0;
 elseif strcmp(ext,'m') | isempty(ext),
	cmd1=['rotfe ' FILE ];	runQ=1;	% run script file
	cmd=['load -mat ' FILE '.mtx'];
 elseif strcmp(ext,'mtx')
	cmd=['load -mat ' FILE '.mtx']; runQ=0;
 elseif ext==[],                                % no extension, assuming template file
	cmd1=['rotfe ' FILE ];	runQ=1;	% run script file
	cmd=['load -mat ' FILE '.mtx'];
 end
 if runQ, eval(cmd1,FAIL); end
 eval(cmd,FAIL); % run script to create model 
end 

% note, the value of Cd is not determined yet,
 if TYPE==0,
  [phi l1 lambda]=modes6(M,Ws*G+D,K,[],Ws==0); 
  end
  if nargout>0, phio=phi; end
    fb=whirldir(phi);
%  perform further operations, i.e. save plot etc. 

 if Save, % save requested
	[dd jj ]=sort(imag(lambda));
	lambda=lambda(jj); phi=phi(:,jj); 
%	if length(l1)>0, l1=l1(:,jj); end
   FILE1=[nameext(FILE) '.mod'];
   FAIL=['error([''rotmodes:  error attempting save, file >> '' FILE1 '' << ''])'];
   cmd=['save ' FILE1 ' phi lambda l1 '];
    eval(cmd,FAIL);
 end
 

  
if Plot, % plot modes
  [name ext]=nameext(FILE);
   rramodes(name);
end 
 
