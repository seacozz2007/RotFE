function Rout=disc(cmd,arg1,arg2)
%shaft
%
% a simple graphical tool for creating complex shafts
% this is  a utility for rotfe
%
% I. bucher 30.8.2000
%
% This program is not in the public doain
% All rights belong to the author Dr. I. Bucher bucher@technion.ac.il

global ROT

if nargin<1
   cmd='start';
end
switch cmd
case 'start'
   shaft
   disc_gui
   set(findobj('tag','shaft_gui'),'vis','off');
   set(findobj('tag','disc_gui'),'vis','on');

   ROT.NODES=[];
   ROT.DISCS=[];
   ROT.ELEMENTS=[];
   ROT.MATERIALS=[];
   axis off
   hold on
case 'set'
   fig=findobj('tag','shaft_gui');
   if isempty(fig), shaft_gui, end
   ROT=arg1;
   draw_elements
case 'delete'
   [E nu rho L do di no z1 z2]=get_fields;
   if nargin>1,   no=arg1;end
   
   ne=size(ROT.ELEMENTS,1);
   if  (no>0) & (ne>no)
      ROT.ELEMENTS(no:end,1:2)=ROT.ELEMENTS(no:end,1:2)-L;
      ROT.ELEMENTS(no,:)=[];
      
      draw_elements;
   end
   
    
case 'draw'
   if nargin<2,    draw_elements;
   else
      draw_elements(arg1);
   end
   
case 'accept'
   [E nu rho  do w di no z ]=get_fields;
   [h_E h_nu h_rho h_do h_w h_di h_no h_z ]=get_handles;
   
   ne=size(ROT.DISCS,1);
   if (no>ne) | (ne==0)
      ROT.DISCS=[   ROT.DISCS;  no do w di E nu rho z];
   else
      ROT.DISCS(no,:)=[no do w di E nu rho z];
   end
   
   %   ROT.NODES=unique(sort([ROT.NODES z1 z2]));
   
   %draw_elements;
   draw_discs
   get_set(   h_z,z);
   get_set(   h_no,no+1);
   
case 'element'
   [h_E h_nu h_rho h_do h_w h_di h_no h_z1 h_z2]=get_handles;
   no=get_set(h_no,arg1);
   if no>0,
      title(sprintf('pressed element # %d',no));  
      disp_element(no);      
      draw_elements(no);
   end
   
case 'save'
   % generate ELEMENTS, MATERIALS, NODES
   NODES=unique(sort([ROT.ELEMENTS(:,1);ROT.ELEMENTS(:,2)]));
   ne=size(ROT.ELEMENTS,1);
   ELEMENTS=[];
   % MATERIALS
   %[	z1 z2 do di E nu rho];
   MAT=unique(ROT.ELEMENTS(:,[5:7]),'rows');
   nm=size(MAT,1);
   for q=1:ne
      n1=ROT.ELEMENTS(q,1);       n2=ROT.ELEMENTS(q,2);
      n1a=find(n1==NODES);       n2a=find(n2==NODES); 
      matno=find(sum(abs(MAT-ones(nm,1)*ROT.ELEMENTS(q,5:7))')==0);
      ELEMENTS=[ELEMENTS;n1a n2a  ROT.ELEMENTS(q,3:4) matno];
   end
   R.MATERIALS=MAT(:,[1 3 2]);
   R.ELEMENTS=ELEMENTS;
   R.NODES=NODES;
   R.DISCS=[];
   R.SPRINGS=[];
   R.BCNodeDir=[];
   R.BCNodeDir=[];
   if nargout>0
      Rout=R;
   else
      assignin('base','R',R);   
      disp(' ')
      disp(' A structure    -->  R')
      disp(' ')
      disp(' containing the Rotfe model geometry')
      disp(' was created in the base workspace of Matlab')
   end
end

function draw_discs(no)
global ROT
kids=get(gca,'child'); delete(kids);
ne=size(ROT.DISCS,1); c='b';
% ROT.DISCS(no,:)=[no do w di E nu rho z];

for q=1:ne
   no=ROT.DISCS(q,1);       w=ROT.DISCS(q,3);
   do=ROT.DISCS(q,2);       di=ROT.DISCS(q,4);
	z=ROT.DISCS(q,8);        
      
   if nargin>0
      if q==no, c='y'; else, c='b'; end
   end
   
   h(1)=fill([z-w/2 z+w/2 z+w/2 z-w/2],[di di do do]/2,c);
   h(1)=fill([z-w/2 z+w/2 z+w/2 z-w/2],-[di di do do]/2,c);
   set(h,'ButtonDownFcn',['disc(''disc''' ',' int2str(q) ')'])
end


function draw_elements(no)
global ROT
kids=get(gca,'child'); delete(kids);
%   [E nu rho L do w di no z ]=get_fields;

ne=size(ROT.ELEMENTS,1); c='r';
for q=1:ne
   z=ROT.ELEMENTS(q,1);       z2=ROT.ELEMENTS(q,2);
   do=ROT.ELEMENTS(q,3);       di=ROT.ELEMENTS(q,4);
   if nargin>0
      if q==no, c='y'; else, c='r'; end
   end
   
   h(1)=fill([z1 z2 z2 z1],[di di do do]/2,c);
   h(2)=fill([z1 z2 z2 z1],-[di di do do]/2,c);
   set(h,'ButtonDownFcn',['disc(''element''' ',' int2str(q) ')'])
end

function disp_element(no)
global ROT
%S={'disc_E' 'disc_nu' 'disc_rho' 'disc_do' 'disc_di' 'disc_no'};
[E nu rho L do di no]=get_fields;
[h_E h_nu h_rho h_L h_do h_di h_no]=get_handles;
get_set(   h_z1,ROT.ELEMENTS(no,1));
get_set(   h_z2,ROT.ELEMENTS(no,2));
get_set(   h_do,ROT.ELEMENTS(no,3));
get_set(   h_di,ROT.ELEMENTS(no,4));
get_set(   h_E,ROT.ELEMENTS(no,5));
get_set(   h_nu,ROT.ELEMENTS(no,6));
get_set(   h_rho,ROT.ELEMENTS(no,7));
get_set(   h_L,z2-z1);

function x=get_set(h,y)
%
try
   s=get(h,'string'); 
   if nargin>1, x=y; else,  x=eval(s); end
   set(h,'string',num2str(x));
catch
   x=-1;
      set(h,'string','err');
end


function varargout=get_fields()
%[h_E,h_nu,h_rho,h_L,h_do,h_di,h_no,h_z1,h_z2]=get_fields()
%
%
S={'disc_E' 'disc_nu' 'disc_rho' 'disc_do' 'disc_w' 'disc_di' 'disc_no' 'disc_z'};
%   [E nu rho  do w di no z ]=get_fields;

nout=nargout;
for q=1:nout
   h=findobj('tag',S{q});
   try
      varargout(q)={eval(get(h,'string'))};
   catch
      varargout(q)={-9999};
   end
end

function varargout=get_handles()
%[h_E,h_nu,h_rho,h_L,h_do,h_di,h_no,h_z1,h_z2]=get_fields()
%
%
S={'disc_E' 'disc_nu' 'disc_rho' 'disc_do' 'disc_w' 'disc_di' 'disc_no'  'disc_z'};
nout=nargout;
for q=1:nout
   h=findobj('tag',S{q});
   try
      varargout(q)={h};
   catch
      varargout(q)={-1};
   end
end

