function [phio,lambda,fb]=tmodeplot(FILE)  
% [phi,lambda,fb]=tmodeplot(FILE) 
%   
%
% INPUTS:
%	FILE:   input file, can be template.m 
%
% EXAMPLES:
%	 [phi,lambda,fb]=rotmodes('template');  % runs rotfe('template'), 
%						% assumes W=0

% By  I. Bucher  8-9-1998
% Rev. 2.0

if nargin<1, FILE='simple1s'; end
tool=mfilename;

if strcmp(FILE,'go'),
   
   fig=findobj('tag',[tool 'fig']);
   hfile=findobj('tag',[tool  'file']);  file=get(hfile,'string');
   Rot=rottors(file);
   hed=findobj('tag',[tool  'modeN']); 
   set(hed,'visible','on');
   hN=findobj('tag',[tool  'N']);  set(hN,'visible','on');
   hcomp=findobj('tag',[tool  'compute']); 
   CoL=get(hcomp,'value');
   [phi lambda]=mkeig(Rot.K,Rot.M); 
   nv=size(phi,1)+length(Rot.TorsionBC);
   phi1=zeros(nv,size(phi,2)); ind=1:nv; ind(Rot.TorsionBC)=[];
   phi1(ind,:)=phi; phi=phi1;
   Rot.eigenvector=phi; Rot.eigenvalue=sqrt(lambda);  
   set(fig,'userdata',Rot);
   Nm=length(lambda);
   S=['1'];  for q=2:Nm, S=str2mat(S,int2str(q)); end
   set(hed,'string',S,'value',1);
   
   tmodeplot('display')
   % <><><><><><><><><><><><>
elseif  strcmp(FILE,'display'),
   hed=findobj('tag',[tool  'modeN']);  n=get(hed,'value');
    
   fig=findobj('tag',[tool 'fig']);
   Rot=get(fig,'userdata');
   ax=findobj('tag',[tool 'ax']);  % get the right axis
   lambda=abs(Rot.eigenvalue);
   drawrot(Rot.File), hold on
   y=Rot.eigenvector(:,n);
   ax1=axis; my=max(abs(ax1(3:4))); y=y*my/max(abs(y));
   plot(Rot.NODES,y,'linewidth',3);
       
   fr=(lambda(n))/2/pi;fr=['f=' num2str(fr) ' Hz'];
   title(['mode ' num2str(n) '   ' fr ]); 
   
   % <><><><>
else, % init
   
    
   fig=findobj('tag',[tool 'fig']);
   if length(fig)==0, 
      fig=figure('unit','normal','pos',...
         [0.0100    0.1    0.6000    0.8],...
         'tag',[tool 'fig'],'menubar','none',...
         'name','Torsional mode shapes',...
         'numbertitle','off'); 
      
   else, figure(gcf), return
   end
   ax=axes; set(ax,'tag',[tool 'ax'])
   set(gca,'pos',[.1 .45 .8 .5]); axis off
    feval('drawrot',FILE);
   % ----------------------------------------------------------------------
   goahead=[tool '(''go'');'];
   uicontrol('style','frame',...
      'foregroundcolor','w',...
      'unit','normal','pos',[.01 .01 .89 .39],...
      'backgroundcolor','b');
   uicontrol('style','text',...
      'string','  Template file (model)',...
      'unit','normal','pos',[.03 .19 .5 .08],'units','normalized',...
      'backgroundcolor','y');
   
   cb=[101 118 97 108 40 91 39 70 73 76 69 61 103 101 116 40 102 105 110 100 111 98 106 40 103 99 102 44 39 39 116 97 103 39 39 44 39 39 116 109 111 100 101 112 108 111 116 102 105 108 101 39 39 41 44 39 39 115 116 114 105 110 103 39 39 41 59 32 100 114 97 119 114 111 116 40 70 73 76 69 41 59 32 97 120 105 115 32 97 117 116 111 39 93 44 91 39 32 116 105 116 108 101 40 39 39 85 115 101 32 101 120 105 115 116 105 110 103 32 84 101 109 112 108 97 116 101 32 70 105 108 101 39 39 41 39 93 41 59]; 
      uicontrol('style','edit',...
      'string',FILE,'units','normalized',...
      'pos',[.55 .19 .15 .08],...
      'tag',[tool 'file'],...
      'callback',setstr(cb)  ,...
      'backgroundcolor','w'); 
     uicontrol('style','text',...
      'string',' show mode number',...
      'units','normalized','pos',[.03 .09 .5 .08],...
      'backgroundcolor','y');
   h1= uicontrol('style','popup',...
      'string','0','units','normalized',...
      'pos',[.55 .11 .15 .08],...
      'tag',[tool 'modeN'],...
      'callback',[tool '(''display'');'],...
      'backgroundcolor','w');
    uicontrol('style','push',...
      'string','GO','units','normalized',...
      'pos',[.72 .03 .16 .1],...
      'callback',goahead);
    
   set([h1 ],'visible','off');
   % ---------------------------------------------------------------------------------
   
end


